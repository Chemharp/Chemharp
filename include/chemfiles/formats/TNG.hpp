// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_FORMAT_TNG_HPP
#define CHEMFILES_FORMAT_TNG_HPP

#include <cstddef>
#include <cstdint>

#include <string>
#include <vector>

#include "chemfiles/File.hpp"
#include "chemfiles/Format.hpp"

#include "chemfiles/files/TNGFile.hpp"

namespace chemfiles {
class Frame;
class FormatMetadata;

/// TNG file format reader.
class TNGFormat final: public Format {
public:
    TNGFormat(std::string path, File::Mode mode, File::Compression compression);

    void read_step(size_t step, Frame& frame) override;
    void read(Frame& frame) override;
    size_t nsteps() override;
private:
    void read_positions(Frame& frame);
    void read_velocities(Frame& frame);
    void read_cell(Frame& frame);
    void read_topology(Frame& frame);

    /// Associated TNG file
    TNGFile tng_;
    /// Scale factor for all length-dependent data:
    /// positions, velocities, forces, and box shape.
    double distance_scale_factor_ = -1;
    /// The next step to read -- in chemfiles numbering
    /// Chemfiles frames are numbered successively: 0, 1, 2, 3; without
    /// accounting for the underlying simulation step
    size_t step_ = 0;
    /// The list of steps in the file -- in TNG numbering
    /// TNG frames are numbered by the MD step: 0, 10, 20, 30
    std::vector<int64_t> tng_steps_;
    /// The number of atoms in the current frame
    int64_t natoms_ = 0;
};

template<> const FormatMetadata& format_metadata<TNGFormat>();

} // namespace chemfiles

#endif
