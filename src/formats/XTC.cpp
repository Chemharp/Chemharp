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
#include "chemfiles/UnitCell.hpp"

#include "chemfiles/files/XDRFile.hpp"
#include "chemfiles/formats/XTC.hpp"

// definitions from the xdrfile library
#define XTC_MAGIC 1995

/* XTC small header size (natoms<=9).
 *  > int(4) magic
 *  > int(4) natoms
 *  > int(4) step
 *  > float(4) time
 *  > 9xfloat(4) box
 *  > int(4) natoms (again)
 */
#define XTC_SMALL_HEADER_SIZE 56

/* Size of uncompressed coordinates for one atom.
 * 3xfloat(4) x
 */
#define XTC_SMALL_COORDS_SIZE 12

/* XTC header size (natoms>=10).
 * Compressed trajectories contain some additional values:
 *  > float(4) precision
 *  > 3xint(4) minint
 *  > 3xint(4) maxint
 *  > int(4) smallidx
 * See `xdrfile_compress_coord_double()`.
 */
#define XTC_HEADER_SIZE (XTC_SMALL_HEADER_SIZE + 32)

using namespace chemfiles;

template <> const FormatMetadata& chemfiles::format_metadata<XTCFormat>() {
    static FormatMetadata metadata;
    metadata.name = "XTC";
    metadata.extension = ".xtc";
    metadata.description = "GROMACS XTC binary format";
    metadata.reference = "http://manual.gromacs.org/current/reference-manual/file-formats.html#xtc";

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

static void get_cell(std::vector<float>& box, const Frame& frame);
static void get_positions(std::vector<float>& x, const Frame& frame);

XTCFormat::XTCFormat(std::string path, File::Mode mode, File::Compression compression)
    : file_(std::move(path), mode) {
    if (compression != File::DEFAULT) {
        throw format_error("XTC format does not support compression");
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

size_t XTCFormat::nsteps() { return frame_offsets_.size(); }

void XTCFormat::read_step(size_t step, Frame& frame) {
    step_ = step;
    file_.seek(frame_offsets_[step_]);
    read(frame);
}

void XTCFormat::read(Frame& frame) {
    FrameHeader header = read_frame_header();

    frame.set_step(header.step);                         // actual step of MD Simulation
    frame.set("time", static_cast<double>(header.time)); // time in pico seconds
    frame.resize(header.natoms);

    const auto box = file_.read_gmx_box();
    frame.set_cell(box);

    size_t natoms_again = file_.read_single_size_as_i32();
    if (natoms_again != header.natoms) {
        throw format_error("contradictory number of atoms in XTC file at '{}': expected {}, got {}",
                           file_.path(), header.natoms, natoms_again);
    }

    std::vector<float> x(header.natoms * 3);
    if (header.natoms <= 9) {
        file_.read_f32(x);
    } else {
        float precision = file_.read_gmx_compressed_floats(x);
        frame.set("xtc_precision", static_cast<double>(precision));
    }
    auto positions = frame.positions();
    assert(x.size() == 3 * positions.size());
    for (size_t i = 0; i < frame.size(); i++) {
        // Factor 10 because the cell lengths are in nm in the XTC format
        positions[i][0] = static_cast<double>(x[i * 3]) * 10.0;
        positions[i][1] = static_cast<double>(x[i * 3 + 1]) * 10.0;
        positions[i][2] = static_cast<double>(x[i * 3 + 2]) * 10.0;
    }

    step_++;
}

XTCFormat::FrameHeader XTCFormat::read_frame_header() {
    try {
        const int32_t magic = file_.read_single_i32();
        if (magic != XTC_MAGIC) {
            throw format_error("invalid XTC file at '{}': "
                               "expected XTC_MAGIC '{}', got '{}'",
                               file_.path(), XTC_MAGIC, magic);
        }

        FrameHeader header = {
            file_.read_single_size_as_i32(), // natoms
            file_.read_single_size_as_i32(), // step
            file_.read_single_f32(),         // time
        };

        return header;
    } catch (const Error& e) {
        throw format_error("could not read XTC header from '{}': {}", file_.path(), e.what());
    }
}

static int32_t round_to_int_boundary(int32_t x) {
    // Rounding to the next 32-bit boundary
    return (x + 3) & ~0x03;
}

void XTCFormat::determine_frame_offsets() {
    uint64_t cur_pos = file_.tell();
    file_.seek(0L);
    FrameHeader header = read_frame_header();

    natoms_ = header.natoms;

    const uint64_t filesize = file_.file_size();

    frame_offsets_.clear();
    frame_offsets_.emplace_back(0);

    // GROMACS does not bother with compression for nine atoms or less
    if (header.natoms <= 9) {
        const uint64_t framebytes =
            static_cast<uint64_t>(XTC_SMALL_HEADER_SIZE + header.natoms * XTC_SMALL_COORDS_SIZE);
        file_.seek(framebytes);

        const uint64_t nframes = filesize / framebytes;
        frame_offsets_.reserve(static_cast<size_t>(nframes));

        for (uint64_t i = 1; i < nframes; ++i) {
            frame_offsets_.emplace_back(i * framebytes);
        }
    } else {
        file_.seek(XTC_HEADER_SIZE);

        uint64_t framebytes = static_cast<uint64_t>(round_to_int_boundary(file_.read_single_i32()));

        const size_t est_nframes = static_cast<size_t>(filesize / (framebytes + XTC_HEADER_SIZE));
        frame_offsets_.reserve(est_nframes);

        while (true) {
            file_.skip(framebytes + XTC_HEADER_SIZE);

            const uint64_t frame_pos = file_.tell() - XTC_HEADER_SIZE;
            try {
                framebytes = static_cast<uint64_t>(round_to_int_boundary(file_.read_single_i32()));
            } catch (const Error&) {
                break;
            }
            frame_offsets_.emplace_back(frame_pos);
        }
    }

    file_.seek(cur_pos);
}

void XTCFormat::write(const Frame& frame) {
    const size_t natoms = frame.size();
    if (frame_offsets_.empty() && step_ == 0) {
        natoms_ = natoms;
    } else if (natoms_ != natoms) {
        throw format_error(
            "XTC format does not support varying numbers of atoms: expected {}, but got {}",
            natoms_, natoms);
    }

    FrameHeader header = {
        natoms,                                                          // natoms
        frame.step(),                                                    // step
        static_cast<float>(frame.get("time").value_or(0.0).as_double()), // time
    };
    write_frame_header(header);

    std::vector<float> box(9);
    get_cell(box, frame);
    file_.write_f32(box);

    file_.write_single_i32(static_cast<int32_t>(header.natoms)); // natoms (again)

    std::vector<float> x(natoms * 3);
    get_positions(x, frame);
    if (natoms <= 9) {
        file_.write_f32(x);
    } else {
        const float precision =
            static_cast<float>(frame.get("xtc_precision").value_or(1000.0).as_double());
        file_.write_gmx_compressed_floats(x, precision);
    }

    step_++;
}

void XTCFormat::write_frame_header(const FrameHeader& header) {
    file_.write_single_i32(XTC_MAGIC);

    file_.write_single_i32(static_cast<int32_t>(header.natoms));
    file_.write_single_i32(static_cast<int32_t>(header.step));
    file_.write_single_f32(static_cast<float>(header.time));
}

void get_cell(std::vector<float>& box, const Frame& frame) {
    assert(box.size() == 9);
    // Factor 10 because the lengths are in nm in the XTC format
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
        // Factor 10 because the lengths are in nm in the XTC format
        x[i * 3] = static_cast<float>(positions[i][0] / 10.0);
        x[i * 3 + 1] = static_cast<float>(positions[i][1] / 10.0);
        x[i * 3 + 2] = static_cast<float>(positions[i][2] / 10.0);
    }
}
