// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_TNG_FORMAT_HPP
#define CHEMFILES_TNG_FORMAT_HPP

#include "chemfiles/Format.hpp"
#include "chemfiles/files/TNGFile.hpp"

namespace chemfiles {

/// [TNG][TNG] file format reader.
///
/// [TNG]: http://dx.doi.org/10.1007/s00894-010-0948-5
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

    /// Reference to the associated file
    TNGFile tng_;
    /// Scale factor for all lenght dependent data e.g. positions, velocities, forces, box shape
    float distance_scale_factor_ = -1;
    /// The next step to read
    int64_t step_ = 0;
    /// The number of atoms in the current frame
    int64_t natoms_ = 0;
};

template<> FormatInfo format_information<TNGFormat>();

} // namespace chemfiles

#endif
