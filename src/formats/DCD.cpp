// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>

#include <string>
#include <utility>
#include <memory>
#include <set>
#include <vector>

#include "chemfiles/File.hpp"
#include "chemfiles/Property.hpp"
#include "chemfiles/types.hpp"
#include "chemfiles/warnings.hpp"
#include "chemfiles/error_fmt.hpp"

#include "chemfiles/Frame.hpp"
#include "chemfiles/UnitCell.hpp"
#include "chemfiles/FormatMetadata.hpp"

#include "chemfiles/files/BinaryFile.hpp"
#include "chemfiles/formats/DCD.hpp"

static double cos_to_angle_degrees(double cos) {
    constexpr double pi = 3.141592653589793238463;
    // use arcsin instead of acos to ensure a cos of 0 is
    // properly mapped to 90 degrees
    return 90.0 - std::asin(cos) * 90.0 / (pi / 2.0);
}

static size_t checked_cast(int32_t value) {
    if (value < 0) {
        throw chemfiles::format_error(
            "invalid value in DCD file: expected a positive integer, got {}",
            value
        );
    }
    return static_cast<size_t>(value);
}

static size_t checked_cast(int64_t value) {
    if (value < 0) {
        throw chemfiles::format_error(
            "invalid value in DCD file: expected a positive integer, got {}",
            value
        );
    }
    return static_cast<size_t>(value);
}


using namespace chemfiles;

static std::unique_ptr<BinaryFile> open_dcd_file(std::string path, File::Mode mode, bool& use_64_bit_markers) {
    if (mode == File::WRITE) {
        return BinaryFile::open_native(std::move(path), mode);
    }

    auto file = LittleEndianFile(path, mode);

    if (mode == File::APPEND && file.file_size() == 0) {
        return BinaryFile::open_native(std::move(path), mode);
    }

    // we need to check multiple variants of the DCD format: 32 and 64-bits
    // fortran record markers, and little or big endian files.
    //
    // The file should start with 84, and then the magic 'CORD'.
    file.seek(0);

    char data[8] = {0};
    file.read_char(data, 8);

    if (data[0] == 84 && data[1] == 0 && data[2] == 0 && data[3] == 0) {
        if (data[4] == 'C' && data[5] == 'O' && data[6] == 'R' && data[7] == 'D') {
            use_64_bit_markers = false;
            return std::make_unique<LittleEndianFile>(std::move(file));
        } else if (data[4] == 0 && data[5] == 0 && data[6] == 0 && data[7] == 0) {
            // We might be using 64-bit record markers, check for CORD
            char extra[4] = {0};
            file.read_char(extra, 4);
            if (extra[0] == 'C' && extra[1] == 'O' && extra[2] == 'R' && extra[3] == 'D') {
                use_64_bit_markers = true;
                return std::make_unique<LittleEndianFile>(std::move(file));
            }
        }
    } else if (data[0] == 0 && data[1] == 0 && data[2] == 0) {
        if (data[3] == 84 && data[4] == 'C' && data[5] == 'O' && data[6] == 'R' && data[7] == 'D') {
            use_64_bit_markers = false;
            return std::make_unique<BigEndianFile>(path, mode);
        } else if (data[3] == 0 && data[4] == 0 && data[5] == 0 && data[6] == 0 && data[7] == 84) {
            // We might be using 64-bit record markers, check for CORD
            char extra[4] = {0};
            file.read_char(extra, 4);
            if (extra[0] == 'C' && extra[1] == 'O' && extra[2] == 'R' && extra[3] == 'D') {
                use_64_bit_markers = true;
                return std::make_unique<BigEndianFile>(path, mode);
            }
        }
    }

    // if we reach here, the file is not a DCD file
    throw format_error(
        "this is not a DCD file, first 8 bytes are 0x{:x}{:x}{:x}{:x} 0x{:x}{:x}{:x}{:x}",
        data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7]
    );
}


DCDFormat::DCDFormat(std::string path, File::Mode mode, File::Compression compression):
    file_(nullptr)
{
    if (compression != File::DEFAULT) {
        throw format_error("compression is not supported for DCD files");
    }

    try {
        file_ = open_dcd_file(path, mode, options_.use_64_bit_markers);
    } catch (const Error& e) {
        throw format_error("unable to open '{}': {}", path, e.what());
    }

    if (mode == File::WRITE || (mode == File::APPEND && file_->file_size() == 0)) {
        return;
    }

    file_->seek(0);
    this->read_header();

    if (mode == File::APPEND) {
        if (options_.has_4d_data) {
            throw format_error("can not append to file with 4D data");
        }

        if (options_.use_64_bit_markers) {
            throw format_error("can not append to file with 64-bit markers");
        }

        if (!fixed_atoms_.empty()) {
            throw format_error("can not append to file with fixed atoms");
        }

        file_->seek(file_->file_size());
    }
}

size_t DCDFormat::nsteps() {
    return n_frames_;
}

void DCDFormat::read(Frame& frame) {
    this->read_step(step_, frame);
    step_++;
}

void DCDFormat::read_step(size_t step, Frame& frame) {
    step_ = step;

    if (step_ == 0) {
        file_->seek(header_size_);
    } else {
        file_->seek(header_size_ + first_frame_size_ + (step_ - 1) * frame_size_);
    }

    frame.set_cell(this->read_cell());
    read_positions(frame);

    // set frame properties
    if (timesteps_.dt != 0.0 && timesteps_.step != 0) {
        auto simulation_step = static_cast<double>(timesteps_.step * step_ + timesteps_.start);
        frame.set("time", timesteps_.dt * simulation_step);
    }

    if (!title_.empty()) {
        frame.set("title", title_);
    }
}

size_t DCDFormat::read_marker() {
    if (options_.use_64_bit_markers) {
        return checked_cast(file_->read_single_i64());
    } else {
        return checked_cast(file_->read_single_i32());
    }
}

void DCDFormat::expect_marker(size_t size) {
    auto marker = this->read_marker();
    if (marker != size) {
        throw format_error("invalid fortran record marker, expected {} got {}", size, marker);
    }
}

void DCDFormat::read_header() {
    auto header_size = this->read_marker();
    assert(header_size == 84);

    // first, check the format using the last integer in the header
    auto header_start = file_->tell();

    file_->seek(header_start + 80);
    options_.charmm_version = static_cast<uint8_t>(file_->read_single_i32());
    options_.charmm_format = options_.charmm_version != 0;

    // then read the rest of the header (the first 4 bytes are CORD)
    file_->seek(header_start + 4);
    n_frames_ = checked_cast(file_->read_single_i32());
    timesteps_.start = checked_cast(file_->read_single_i32());
    timesteps_.step = checked_cast(file_->read_single_i32());

    // 20 unused bytes in the header
    file_->seek(file_->tell() + 20);
    auto n_fixed_atoms = checked_cast(file_->read_single_i32());

    if (options_.charmm_format) {
        timesteps_.dt = static_cast<double>(file_->read_single_f32());

        if (file_->read_single_i32() != 0) {
            options_.charmm_unitcell = true;
        }

        if (file_->read_single_i32() == 1) {
            options_.has_4d_data = true;
        }
    } else {
        // `timesteps_.dt` uses a 64-bite float in X-PLOR format
        timesteps_.dt = file_->read_single_f64();
    }

    file_->seek(header_start + 84);
    this->expect_marker(84);

    auto title_size = this->read_marker();
    if (title_size < 4 || (title_size - 4) % 80 != 0) {
        if (title_size != 0) {
            warning("DCD reader", "invalid title record size ({}), skipping title section", title_size);
            file_->seek(file_->tell() + static_cast<uint64_t>(title_size));
        }
    } else {
        auto n_lines = checked_cast(file_->read_single_i32());
        if (n_lines != (title_size - 4) / 80) {
            warning(
                "DCD reader",
                "title record size ({}) does not match the number of title lines ({}), skipping title section",
                title_size, n_lines
            );
            file_->seek(file_->tell() + static_cast<uint64_t>(title_size - 4));
        }

        auto buffer = std::vector<char>(title_size - 4, '\0');
        file_->read_char(buffer);

        // each line in the title might be padded with NULL or space characters,
        // or just NULL-terminated and contain garbage. This code tries to clean
        // up the title, removing the risk of reading in garbage.
        for (size_t line_i=0; line_i<n_lines; line_i++) {
            for (size_t char_i=line_i * 80; char_i<(line_i + 1) * 80; char_i++) {
                if (buffer[char_i] == '\0') {
                    title_ += '\n';
                    break;
                } else {
                    title_ += buffer[char_i];
                }
            }
        }
    }
    this->expect_marker(title_size);

    // the next record contains the number of atoms
    this->expect_marker(sizeof(int32_t));
    n_atoms_ = checked_cast(file_->read_single_i32());
    buffer_.resize(n_atoms_);
    this->expect_marker(sizeof(int32_t));

    n_free_atoms_ = n_atoms_;
    if (n_fixed_atoms != 0) {
        n_free_atoms_ = n_atoms_ - n_fixed_atoms;
        this->expect_marker(sizeof(int32_t) * n_free_atoms_);

        auto free_atoms = std::vector<int32_t>(n_free_atoms_);
        file_->read_i32(free_atoms);

        this->expect_marker(sizeof(int32_t) * n_free_atoms_);

        // translate from fortran 1-based indexes to 0-based indexes for
        // free atoms
        std::for_each(std::begin(free_atoms), std::end(free_atoms), [&](int32_t& index) {
            assert(index >= 1 && index <= static_cast<int32_t>(n_atoms_));
            index -= 1;
        });

        // make sure the free atoms are sorted before doing a binary search
        // with `std::lower_bound` below
        std::sort(std::begin(free_atoms), std::end(free_atoms));

        fixed_atoms_.reserve(n_atoms_);
        for (int32_t i=0; i<static_cast<int32_t>(n_atoms_); i++) {
            auto data = fixed_atom_data_t{};
            auto it = std::lower_bound(std::begin(free_atoms), std::end(free_atoms), i);
            if (it != std::end(free_atoms) && *it == i) {
                // this is a free atom
                data.fixed = false;
                data.free_index = static_cast<size_t>(it - std::begin(free_atoms));
                fixed_atoms_.push_back(data);
            } else {
                // this is a fixed atom
                data.fixed = true;
                data.fixed_coord = Vector3D();
                fixed_atoms_.push_back(data);
            }
        }
    }

    // we are done with the header, now we can compute the size of a frame and
    // use it to check if the file size and number of frames in the header match
    header_size_ = file_->tell();

    auto compute_frame_size = [&](size_t n_atoms) {
        size_t result = 0;

        size_t marker_size = options_.use_64_bit_markers ? 8 : 4;
        if (options_.charmm_format && options_.charmm_unitcell) {
            // 6 doubles and two markers for the unit cell
            result += 2 * marker_size + 6 * sizeof(double);
        }
        // three coordinates sets, each with two markers
        result += 3 * (2 * marker_size + n_atoms * sizeof(float));
        if (options_.has_4d_data) {
            // optional fourth dimension
            result += 2 * marker_size + n_atoms * sizeof(float);
        }

        return result;
    };

    first_frame_size_ = compute_frame_size(n_atoms_);
    if (fixed_atoms_.empty()) {
        frame_size_ = first_frame_size_;
    } else {
        frame_size_ = compute_frame_size(n_free_atoms_);
    }

    auto file_size = file_->file_size();
    auto n_frames_from_size = (file_size - header_size_ - first_frame_size_) / frame_size_ + 1;
    if (n_frames_from_size != n_frames_) {
        warning(
            "DCD reader",
            "the file header claims {} frames, but the file size indicate we have {}",
            n_frames_,
            n_frames_from_size
        );
        n_frames_ = static_cast<size_t>(n_frames_from_size);
    }

    if (!fixed_atoms_.empty()) {
        this->read_fixed_coordinates();
    }
}

UnitCell DCDFormat::read_cell() {
    if (!options_.charmm_format || !options_.charmm_unitcell) {
        return {};
    }

    this->expect_marker(6 * sizeof(double));
    double buffer[6] = {0};
    file_->read_f64(buffer, 6);
    this->expect_marker(6 * sizeof(double));

    if (options_.charmm_format && options_.charmm_version > 25) {
        // recent charmm version store cell vectors directly. It seems to be
        // doing this since CHARMM 22, but the readdcd code pretends to be
        // CHARMM 24 most of the time, so we need to check for a higher version.
        auto matrix = Matrix3D(
            buffer[0], buffer[1], buffer[3],
            buffer[1], buffer[2], buffer[4],
            buffer[3], buffer[4], buffer[5]
        );
        return UnitCell(matrix);
    }

    auto lengths = Vector3D(buffer[0], buffer[2], buffer[5]);
    auto alpha = buffer[4];
    auto beta = buffer[3];
    auto gamma = buffer[1];

    if (std::fabs(alpha) <= 1.0 && std::fabs(beta) <= 1.0 && std::fabs(gamma) <= 1.0) {
        // cell angles can be stored either directly or as the cos of the angle.
        // if all of the angles are smaller than 1, we assume they are stored as
        // their cosine
        alpha = cos_to_angle_degrees(alpha);
        beta = cos_to_angle_degrees(beta);
        gamma = cos_to_angle_degrees(gamma);
    }

    auto angles = Vector3D(alpha, beta, gamma);
    return UnitCell(lengths, angles);
}

void DCDFormat::read_positions(Frame& frame) {
    frame.resize(n_atoms_);
    auto positions = frame.positions();

    auto n_atoms_to_read = n_atoms_;
    if (!fixed_atoms_.empty()) {
        if (step_ != 0) {
            n_atoms_to_read = n_free_atoms_;
            for (size_t i=0; i<frame.size(); i++) {
                if (fixed_atoms_[i].fixed) {
                    positions[i] = fixed_atoms_[i].fixed_coord;
                }
            }
        }
    }

    buffer_.resize(n_atoms_to_read);

    // read the X coordinate
    this->expect_marker(sizeof(float) * n_atoms_to_read);
    file_->read_f32(buffer_);
    this->expect_marker(sizeof(float) * n_atoms_to_read);

    if (n_atoms_to_read == n_atoms_) {
        for (size_t i=0; i<n_atoms_; i++) {
            positions[i][0] = static_cast<double>(buffer_[i]);
        }
    } else {
        for (size_t i=0; i<frame.size(); i++) {
            if (!fixed_atoms_[i].fixed) {
                positions[i][0] = static_cast<double>(buffer_[fixed_atoms_[i].free_index]);
            }
        }
    }

    // read the Y coordinate
    this->expect_marker(sizeof(float) * n_atoms_to_read);
    file_->read_f32(buffer_);
    this->expect_marker(sizeof(float) * n_atoms_to_read);

    if (n_atoms_to_read == n_atoms_) {
        for (size_t i=0; i<n_atoms_; i++) {
            positions[i][1] = static_cast<double>(buffer_[i]);
        }
    } else {
        for (size_t i=0; i<frame.size(); i++) {
            if (!fixed_atoms_[i].fixed) {
                positions[i][1] = static_cast<double>(buffer_[fixed_atoms_[i].free_index]);
            }
        }
    }

    // read the Z coordinate
    this->expect_marker(sizeof(float) * n_atoms_to_read);
    file_->read_f32(buffer_);
    this->expect_marker(sizeof(float) * n_atoms_to_read);

    if (n_atoms_to_read == n_atoms_) {
        for (size_t i=0; i<n_atoms_; i++) {
            positions[i][2] = static_cast<double>(buffer_[i]);
        }
    } else {
        for (size_t i=0; i<frame.size(); i++) {
            if (!fixed_atoms_[i].fixed) {
                positions[i][2] = static_cast<double>(buffer_[fixed_atoms_[i].free_index]);
            }
        }
    }

    if (options_.has_4d_data) {
        // skip the fourth dimension
        this->expect_marker(sizeof(float) * n_atoms_to_read);
        file_->seek(file_->tell() + sizeof(float) * n_atoms_to_read);
        this->expect_marker(sizeof(float) * n_atoms_to_read);
    }
}

void DCDFormat::read_fixed_coordinates() {
    auto frame = Frame();
    this->read_step(0, frame);
    assert(fixed_atoms_.size() == frame.size());

    auto positions = frame.positions();
    for (size_t i=0; i<frame.size(); i++) {
        if (fixed_atoms_[i].fixed) {
            fixed_atoms_[i].fixed_coord = positions[i];
        }
    }
}

/******************************************************************************/

void DCDFormat::write(const Frame& frame) {
    if (n_frames_ == 0) {
        // initialize data that will be constant for this file
        n_atoms_ = frame.size();
        n_free_atoms_ = frame.size();

        if (n_atoms_ == 0) {
            throw file_error("can not write a frame with 0 atoms");
        }

        title_ = frame.get<Property::STRING>("title").value_or("");

        options_.charmm_format = true;
        options_.charmm_version = 24;
        options_.charmm_unitcell = frame.cell().shape() != UnitCell::INFINITE;
        options_.use_64_bit_markers = false;
        options_.has_4d_data = false;

        file_->seek(0);
        this->write_header();

        header_size_ = file_->tell();

        frame_size_ = 0;
        size_t marker_size = 4;
        // 6 doubles and two markers for the unit cell
        frame_size_ += 2 * marker_size + 6 * sizeof(double);
        // three coordinates sets, each with two markers
        frame_size_ += 3 * (2 * marker_size + n_atoms_ * sizeof(float));

        first_frame_size_ = frame_size_;
    } else if (n_atoms_ != frame.size()) {
        throw format_error(
            "this file was initialized with {} atoms, can not write a frame with {} atoms to it",
            n_atoms_, frame.size()
        );
    } else if (n_atoms_ != n_free_atoms_) {
        throw format_error("can not append to a file with fixed atoms");
    } else if (options_.has_4d_data) {
        throw format_error("can not append to a file with 4D data");
    } else if (options_.use_64_bit_markers) {
        throw format_error("can not append to a file with 64 bit markers");
    }

    auto title = frame.get<Property::STRING>("title");
    if (title && title.value() != title_) {
        warning("DCD writer",
            "the title of this frame doesn't match the title of the file, "
            "the frame title will be ignored"
        );
    }

    this->write_cell(frame.cell());
    this->write_positions(frame);

    n_frames_ += 1;
    step_ += 1;

    // update the number of frames in the header
    auto current = file_->tell();
    file_->seek(8); // the number of frames is always at offset 8 (1 marker + CORD)
    file_->write_single_i32(static_cast<int32_t>(n_frames_));
    file_->seek(current);
}

void DCDFormat::write_header() {
    this->write_marker(84);

    file_->write_char("CORD", 4);
    file_->write_single_i32(static_cast<int32_t>(n_frames_));
    file_->write_single_i32(static_cast<int32_t>(timesteps_.start));
    file_->write_single_i32(static_cast<int32_t>(timesteps_.step));

    // 16 unused bytes
    file_->write_char("\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0", 16);

    file_->write_single_i32(/* n_degree_of_freedom */ 3 * static_cast<int32_t>(n_atoms_));
    file_->write_single_i32(/* n_fixed_atoms */ 0);
    // TODO: this might lose information: when writing to a file, this is always
    // 0, even if the frame has a *time* property. If we read back, the
    // resulting frame will not have a *time* property at all.
    file_->write_single_f32(static_cast<float>(timesteps_.dt));

    file_->write_single_i32(options_.charmm_unitcell ? 1 : 0);
    file_->write_single_i32(/* has_4d_data */ 0);

    // 28 unused bytes
    file_->write_char("\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0", 28);
    file_->write_single_i32(/* charmm version */ 24);
    this->write_marker(84);

    if (title_.empty()) {
        this->write_marker(sizeof(int32_t));
        file_->write_single_i32(0);
        this->write_marker(sizeof(int32_t));
    } else {
        std::string title = title_;
        if (title.size() % 80 != 0) {
            auto padded_size = (title.size() / 80 + 1) * 80;
            title.resize(padded_size, '\0');
        }
        this->write_marker(title.size() + sizeof(int32_t));
        file_->write_single_i32(static_cast<int32_t>(title.size() / 80));
        file_->write_char(title.data(), title.size());
        this->write_marker(title.size() + sizeof(int32_t));
    }

    this->write_marker(sizeof(int32_t));
    file_->write_single_i32(static_cast<int32_t>(n_atoms_));
    this->write_marker(sizeof(int32_t));
}


void DCDFormat::write_cell(const UnitCell& cell) {
    if (cell.shape() == UnitCell::INFINITE) {
        if (options_.charmm_unitcell) {
            warning("DCD writer",
                "this file contains unit cell information, but we have an "
                "infinite cell, we'll write zeros for the cell lengths"
            );
        } else {
            return;
        }
    } else {
        if (!options_.charmm_unitcell) {
            warning("DCD writer",
                "this file does not store unit cell information, "
                "we'll skip writing the cell"
            );
            return;
        }
    }

    if (!chemfiles::private_details::is_upper_triangular(cell.matrix())) {
        warning("DCD writer",
            "the unit cell is not upper-triangular, positions might not align "
            "with the cell in the file"
        );
    }

    auto lengths = cell.lengths();
    auto angles = cell.angles();

    this->write_marker(6 * sizeof(double));
    double buffer[6] = {
        lengths[0],
        angles[2],
        lengths[1],
        angles[1],
        angles[0],
        lengths[2],
    };
    file_->write_f64(buffer, 6);
    this->write_marker(6 * sizeof(double));
}


void DCDFormat::write_positions(const Frame& frame) {
    const auto& positions = frame.positions();

    buffer_.resize(n_atoms_);
    for (size_t i=0; i<n_atoms_; i++) {
        buffer_[i] = static_cast<float>(positions[i][0]);
    }
    this->write_marker(sizeof(float) * n_atoms_);
    file_->write_f32(buffer_);
    this->write_marker(sizeof(float) * n_atoms_);


    for (size_t i=0; i<n_atoms_; i++) {
        buffer_[i] = static_cast<float>(positions[i][1]);
    }
    this->write_marker(sizeof(float) * n_atoms_);
    file_->write_f32(buffer_);
    this->write_marker(sizeof(float) * n_atoms_);


    for (size_t i=0; i<n_atoms_; i++) {
        buffer_[i] = static_cast<float>(positions[i][2]);
    }
    this->write_marker(sizeof(float) * n_atoms_);
    file_->write_f32(buffer_);
    this->write_marker(sizeof(float) * n_atoms_);
}


void DCDFormat::write_marker(size_t size) {
    if (options_.use_64_bit_markers) {
        file_->write_single_i64(static_cast<int64_t>(size));
    } else {
        file_->write_single_i32(static_cast<int32_t>(size));
    }
}

/******************************************************************************/

template <> const FormatMetadata& chemfiles::format_metadata<DCDFormat>() {
    static FormatMetadata metadata;
    metadata.name = "DCD";
    metadata.extension = ".dcd";
    metadata.description = "DCD binary format";
    metadata.reference = "https://web.archive.org/web/20070406065433/http://www.bio.unizh.ch/docu/acc_docs/doc/charmm_principles/Ch04_mol_dyn.FM5.html";

    metadata.read = true;
    metadata.write = true;
    metadata.memory = false;

    metadata.positions = true;
    metadata.velocities = false;
    metadata.unit_cell = true;
    metadata.atoms = false;
    metadata.bonds = false;
    metadata.residues = false;
    return metadata;
}
