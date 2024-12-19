// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cassert>
#include <cstddef>
#include <cstdint>

#include <array>
#include <string>
#include <utility>
#include <vector>

#include "chemfiles/error_fmt.hpp"
#include "chemfiles/external/optional.hpp"
#include "chemfiles/external/span.hpp"
#include "chemfiles/types.hpp"

#include "chemfiles/File.hpp"
#include "chemfiles/FormatMetadata.hpp"
#include "chemfiles/Frame.hpp"
#include "chemfiles/Property.hpp"
#include "chemfiles/UnitCell.hpp"

#include "chemfiles/files/XDRFile.hpp"
#include "chemfiles/formats/TRR.hpp"

// definitions from the xdrfile library
#define TRR_MAGIC 1993
#define TRR_VERSION "GMX_trn_file"

/* Minimum TRR header size.
 *  > int(4) magic
 *  > string version = uint(4) + uint(4) n + bytes(n+1)
 *  > 10xint(4) ir_size, e_size, box_size, vir_size, pres_size,
 *               top_size, sym_size, x_size, v_size, f_size
 *  > int(4) natoms
 *  > int(4) step
 *  > int(4) nre
 *  > float(4)/double(8) t
 *  > float(4)/double(8) lamda
 * For an empty version string (n=0) this adds up to 72 bytes.
 * Default version string is "GMX_trn_file" with n=12, so 84 bytes are typical.
 * It can have 8 bytes more if we have double time and lambda.
 */
#define TRR_MIN_HEADER_SIZE 72

using namespace chemfiles;

template <> const FormatMetadata& chemfiles::format_metadata<TRRFormat>() {
    static FormatMetadata metadata;
    metadata.name = "TRR";
    metadata.extension = ".trr";
    metadata.description = "GROMACS TRR binary format";
    metadata.reference = "http://manual.gromacs.org/current/reference-manual/file-formats.html#trr";

    metadata.read = true;
    metadata.write = true;
    metadata.memory = false;

    metadata.positions = true;
    metadata.velocities = true;
    metadata.unit_cell = true;
    metadata.atoms = false;
    metadata.bonds = false;
    metadata.residues = false;
    return metadata;
}

static void get_cell(std::vector<float>& box, const Frame& frame);
static void get_positions(std::vector<float>& x, const Frame& frame);
static void get_velocities(std::vector<float>& v, const Frame& frame);

TRRFormat::TRRFormat(std::string path, File::Mode mode, File::Compression compression)
    : file_(std::move(path), mode) {
    if (compression != File::DEFAULT) {
        throw format_error("TRR format does not support compression");
    }

    if (mode == File::READ) {
        determine_frame_offsets();
    } else if (mode == File::APPEND) {
        try {
            determine_frame_offsets();
        } catch (const Error&) {
            // Ignore exceptions, because the file might not exist. If it does,
            // we need to get the number of atoms and frames for appending.
        }
    }
}

size_t TRRFormat::nsteps() { return frame_offsets_.size(); }

void TRRFormat::read_step(size_t step, Frame& frame) {
    step_ = step;
    file_.seek(frame_offsets_[step_]);
    read(frame);
}

void TRRFormat::read(Frame& frame) {
    FrameHeader header = read_frame_header();

    bool has_box = (header.box_size > 0);
    bool has_positions = (header.x_size > 0);
    bool has_velocities = (header.v_size > 0);

    frame.set_step(header.step);            // actual step of MD Simulation
    frame.set("time", header.time);         // time in pico seconds
    frame.set("trr_lambda", header.lambda); // coupling parameter for free energy methods
    frame.set("has_positions", has_positions);
    frame.resize(header.natoms);

    if (has_box) {
        const auto box = file_.read_gmx_box(header.use_double);
        frame.set_cell(box);
    }

    size_t legacy_size = header.vir_size + header.pres_size;
    if (legacy_size > 0) {
        file_.skip(static_cast<uint64_t>(legacy_size));
    }

    if (header.use_double) {
        std::vector<double> dx(header.natoms * 3);
        if (has_positions) {
            file_.read_f64(dx);
            auto positions = frame.positions();
            assert(dx.size() == 3 * positions.size());
            for (size_t i = 0; i < frame.size(); i++) {
                // Factor 10 because the cell lengths are in nm in the TRR format
                positions[i][0] = dx[i * 3] * 10.0;
                positions[i][1] = dx[i * 3 + 1] * 10.0;
                positions[i][2] = dx[i * 3 + 2] * 10.0;
            }
        }
        if (has_velocities) {
            file_.read_f64(dx);
            frame.add_velocities();
            auto velocities = *frame.velocities();
            assert(dx.size() == 3 * velocities.size());
            for (size_t i = 0; i < frame.size(); i++) {
                // Factor 10 because the lengths are in nm in the TRR format
                velocities[i][0] = dx[i * 3] * 10.0;
                velocities[i][1] = dx[i * 3 + 1] * 10.0;
                velocities[i][2] = dx[i * 3 + 2] * 10.0;
            }
        }
    } else {
        std::vector<float> dx(header.natoms * 3);
        if (has_positions) {
            file_.read_f32(dx);
            auto positions = frame.positions();
            assert(dx.size() == 3 * positions.size());
            for (size_t i = 0; i < frame.size(); i++) {
                // Factor 10 because the cell lengths are in nm in the TRR format
                positions[i][0] = static_cast<double>(dx[i * 3]) * 10.0;
                positions[i][1] = static_cast<double>(dx[i * 3 + 1]) * 10.0;
                positions[i][2] = static_cast<double>(dx[i * 3 + 2]) * 10.0;
            }
        }
        if (has_velocities) {
            file_.read_f32(dx);
            frame.add_velocities();
            auto velocities = *frame.velocities();
            assert(dx.size() == 3 * velocities.size());
            for (size_t i = 0; i < frame.size(); i++) {
                // Factor 10 because the lengths are in nm in the TRR format
                velocities[i][0] = static_cast<double>(dx[i * 3]) * 10.0;
                velocities[i][1] = static_cast<double>(dx[i * 3 + 1]) * 10.0;
                velocities[i][2] = static_cast<double>(dx[i * 3 + 2]) * 10.0;
            }
        }
    }

    if (header.f_size > 0) {
        file_.skip(static_cast<uint64_t>(header.f_size));
    }

    step_++;
}

TRRFormat::FrameHeader TRRFormat::read_frame_header() {
    try {
        const int32_t magic = file_.read_single_i32();
        if (magic != TRR_MAGIC) {
            throw format_error("invalid TRR file at '{}': "
                               "expected TRR_MAGIC '{}', got '{}'",
                               file_.path(), TRR_MAGIC, magic);
        }

        const std::string version = file_.read_gmx_string();
        if (version != TRR_VERSION) {
            throw format_error("invalid TRR file at '{}': "
                               "expected TRR_VERSION '{}', got '{}'",
                               file_.path(), TRR_VERSION, version);
        }

        FrameHeader header = {
            false,                           // use_double
            file_.read_single_size_as_i32(), // ir_size
            file_.read_single_size_as_i32(), // e_size
            file_.read_single_size_as_i32(), // box_size
            file_.read_single_size_as_i32(), // vir_size
            file_.read_single_size_as_i32(), // pres_size
            file_.read_single_size_as_i32(), // top_size
            file_.read_single_size_as_i32(), // sym_size
            file_.read_single_size_as_i32(), // x_size
            file_.read_single_size_as_i32(), // v_size
            file_.read_single_size_as_i32(), // f_size

            file_.read_single_size_as_i32(), // natoms
            file_.read_single_size_as_i32(), // step
            file_.read_single_size_as_i32(), // nre
            0.0,                             // time
            0.0,                             // lambda
        };

        // determine real representation (float or double)
        size_t nflsize = 0;
        if (header.box_size > 0) {
            nflsize = header.box_size / 9;
        } else if (header.natoms > 0) {
            if (header.x_size > 0) {
                nflsize = header.x_size / (header.natoms * 3);
            } else if (header.v_size > 0) {
                nflsize = header.v_size / (header.natoms * 3);
            } else if (header.f_size > 0) {
                nflsize = header.f_size / (header.natoms * 3);
            }
        } else {
            throw format_error("invalid TRR file at '{}': "
                               "frame with step '{}' has no data",
                               file_.path(), header.step);
        }
        if (nflsize != sizeof(float) && nflsize != sizeof(double)) {
            throw format_error("invalid TRR file at '{}': "
                               "frame with step '{}' has an invalid size",
                               file_.path(), header.step);
        }
        header.use_double = (nflsize == sizeof(double));

        if (header.use_double) {
            header.time = file_.read_single_f64();
            header.lambda = file_.read_single_f64();
        } else {
            header.time = static_cast<double>(file_.read_single_f32());
            header.lambda = static_cast<double>(file_.read_single_f32());
        }

        return header;
    } catch (const Error& e) {
        throw format_error("could not read TRR header from '{}': {}", file_.path(), e.what());
    }
}

void TRRFormat::determine_frame_offsets() {
    uint64_t cur_pos = file_.tell();
    file_.seek(0L);
    FrameHeader header = read_frame_header();

    natoms_ = header.natoms;

    auto calc_framebytes = [&header]() {
        return (
            header.ir_size + header.e_size + header.box_size + header.vir_size + header.pres_size +
            header.top_size + header.sym_size + header.x_size + header.v_size + header.f_size
        );
    };
    uint64_t framebytes = calc_framebytes();

    uint64_t filesize = file_.file_size();
    auto est_nframes = static_cast<size_t>(filesize / (framebytes + TRR_MIN_HEADER_SIZE));

    frame_offsets_.clear();
    frame_offsets_.emplace_back(0);
    frame_offsets_.reserve(est_nframes);

    while (true) {
        file_.skip(framebytes);

        auto frame_pos = file_.tell();
        try {
            header = read_frame_header();
        } catch (const Error&) {
            break;
        }
        frame_offsets_.emplace_back(frame_pos);

        framebytes = calc_framebytes();
    }

    file_.seek(cur_pos);
}

void TRRFormat::write(const Frame& frame) {
    const size_t natoms = frame.size();
    if (frame_offsets_.empty() && step_ == 0) {
        natoms_ = natoms;
    } else if (natoms_ != natoms) {
        throw format_error(
            "TRR format does not support varying numbers of atoms: expected {}, but got {}",
            natoms_, natoms);
    }

    size_t box_size = sizeof(float) * 3 * 3;
    if (frame.cell().shape() == UnitCell::INFINITE) {
        // no box data
        box_size = 0;
    }

    const size_t dx_size = sizeof(float) * natoms * 3;

    size_t x_size = dx_size;
    if (frame.get("has_positions") && !(*frame.get("has_positions")).as_bool()) {
        // no position data
        x_size = 0;
    }

    size_t v_size = dx_size;
    if (!frame.velocities()) {
        // no velocity data
        v_size = 0;
    }

    FrameHeader header = {
        false,    // use_double
        0,        // ir_size
        0,        // e_size
        box_size, // box_size
        0,        // vir_size
        0,        // pres_size
        0,        // top_size
        0,        // sym_size
        x_size,   // x_size
        v_size,   // v_size
        0,        // f_size

        natoms,                                            // natoms
        frame.step(),                                      // step
        0,                                                 // nre
        frame.get("time").value_or(0.0).as_double(),       // time
        frame.get("trr_lambda").value_or(0.0).as_double(), // lambda
    };
    write_frame_header(header);

    std::vector<float> box(9);
    if (box_size > 0) {
        get_cell(box, frame);
        file_.write_f32(box);
    }

    if (x_size > 0 || v_size > 0) {
        std::vector<float> dx(natoms * 3);
        if (x_size > 0) {
            get_positions(dx, frame);
            file_.write_f32(dx);
        }
        if (v_size > 0) {
            get_velocities(dx, frame);
            file_.write_f32(dx);
        }
    }

    step_++;
}

void TRRFormat::write_frame_header(const FrameHeader& header) {
    file_.write_single_i32(TRR_MAGIC);

    file_.write_gmx_string(TRR_VERSION);

    // use_double is not written and has to be inferred when reading
    file_.write_single_i32(static_cast<int32_t>(header.ir_size));
    file_.write_single_i32(static_cast<int32_t>(header.e_size));
    file_.write_single_i32(static_cast<int32_t>(header.box_size));
    file_.write_single_i32(static_cast<int32_t>(header.vir_size));
    file_.write_single_i32(static_cast<int32_t>(header.pres_size));
    file_.write_single_i32(static_cast<int32_t>(header.top_size));
    file_.write_single_i32(static_cast<int32_t>(header.sym_size));
    file_.write_single_i32(static_cast<int32_t>(header.x_size));
    file_.write_single_i32(static_cast<int32_t>(header.v_size));
    file_.write_single_i32(static_cast<int32_t>(header.f_size));

    file_.write_single_i32(static_cast<int32_t>(header.natoms));
    file_.write_single_i32(static_cast<int32_t>(header.step));
    file_.write_single_i32(static_cast<int32_t>(header.nre));
    file_.write_single_f32(static_cast<float>(header.time));
    file_.write_single_f32(static_cast<float>(header.lambda));
}

void get_cell(std::vector<float>& box, const Frame& frame) {
    assert(box.size() == 9);
    // Factor 10 because the lengths are in nm in the TRR format
    auto matrix = frame.cell().matrix() / 10.0;
    box[0] = static_cast<float>(matrix[0][0]);
    box[1] = static_cast<float>(matrix[1][0]);
    box[2] = static_cast<float>(matrix[2][0]);
    box[3] = static_cast<float>(matrix[0][1]);
    box[4] = static_cast<float>(matrix[1][1]);
    box[5] = static_cast<float>(matrix[2][1]);
    box[6] = static_cast<float>(matrix[0][2]);
    box[7] = static_cast<float>(matrix[1][2]);
    box[8] = static_cast<float>(matrix[2][2]);
}

void get_positions(std::vector<float>& x, const Frame& frame) {
    auto positions = frame.positions();
    assert(x.size() == 3 * positions.size());
    for (size_t i = 0; i < frame.size(); ++i) {
        // Factor 10 because the lengths are in nm in the TRR format
        x[i * 3] = static_cast<float>(positions[i][0] / 10.0);
        x[i * 3 + 1] = static_cast<float>(positions[i][1] / 10.0);
        x[i * 3 + 2] = static_cast<float>(positions[i][2] / 10.0);
    }
}

void get_velocities(std::vector<float>& v, const Frame& frame) {
    auto velocities = *frame.velocities();
    assert(v.size() == 3 * velocities.size());
    for (size_t i = 0; i < frame.size(); i++) {
        // Factor 10 because the lengths are in nm in the TRR format
        v[i * 3] = static_cast<float>(velocities[i][0] / 10.0);
        v[i * 3 + 1] = static_cast<float>(velocities[i][1] / 10.0);
        v[i * 3 + 2] = static_cast<float>(velocities[i][2] / 10.0);
    }
}
