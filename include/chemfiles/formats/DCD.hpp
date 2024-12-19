// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_FORMAT_DCD_HPP
#define CHEMFILES_FORMAT_DCD_HPP

#include <cstddef>
#include <cstdint>

#include <memory>
#include <string>
#include <vector>

#include "chemfiles/File.hpp"
#include "chemfiles/Format.hpp"
#include "chemfiles/types.hpp"

#include "chemfiles/files/BinaryFile.hpp"

namespace chemfiles {

class Frame;
class UnitCell;
class Vector3D;
class FormatMetadata;

/// DCD file reader and writer
///
/// The DCD format was originally introduced by CHARMM, and is now used mainly
/// by NAMD. A lot of different analysis software supports it, making it widely
/// used in the bio-molecules community.
///
/// This format is based on fortran's unformatted binary I/O. Each call to
/// `write` in fortran produces a record marker containing the size in bytes of
/// the field, followed by a binary dump of the data, and then the record marker
/// again. Nowadays, record markers are 32-bit most of the time, but gfortran
/// (and maybe others) have an option to use 64-bit markers, and they were the
/// default at some point. The record data and markers are written using the
/// current machine endianness, forcing us to read both little and big endian
/// versions of this format.
///
/// This implementation is able to read a wide variety of variations of this
/// format, including:
///
/// - files created on little & big endian machines, without having to convert
///   the file first to the current machine endianess
/// - files with fixed atoms
/// - files with 4D dynamic
/// - multiple conventions for the unit cell: `a, gamma, b, beta, alpha, c` with
///   the angles stored in degrees or as the cosine of the angle; as well as the
///   full unit cell vectors (this is only used by CHARMM, starting with version
///   22).
///
/// When writing, this format uses a native endianess file, always outputs the
/// unit cell (infinite unit cells being represented with 0), and uses a 3D
/// format without any fixed atoms.
class DCDFormat final: public Format {
public:
    DCDFormat(std::string path, File::Mode mode, File::Compression compression);

    size_t nsteps() override;
    void read(Frame& frame) override;
    void read_step(size_t step, Frame& frame) override;
    void write(const Frame& frame) override;

private:


    /****** low-level function to read fortran unformatted binary files ******/
    // read a single record size marker from the file. Each record (single
    // fortran `write` statement) starts & ends with its size in bytes.
    size_t read_marker();
    // check that the file contains a marker for the expected `size` at the
    // current position
    void expect_marker(size_t size);

    void write_marker(size_t size);

    /************ high level function specialized for DCD format **************/
    /// read & parse the file header
    void read_header();
    UnitCell read_cell();
    void read_positions(Frame& frame);
    void read_fixed_coordinates();

    void write_header();
    void write_cell(const UnitCell& cell);
    void write_positions(const Frame& frame);

    std::unique_ptr<BinaryFile> file_;
    /// which variant of the DCD format are we trying to read?
    struct {
        /// Does the file use 32 or 64-bits fortran record markers?
        bool use_64_bit_markers = false;
        /// Does the file use the CHARMM or X-PLOR variant of the DCD format
        bool charmm_format = false;
        /// The CHARMM version this file says it uses
        uint8_t charmm_version = 0;
        /// Does the file contains unit cell information?
        bool charmm_unitcell = false;
        /// CHARMM supports adding an extra dimension to MD simulation, which is
        /// also saved in DCD files
        bool has_4d_data = false;
    } options_;
    /// end of the header in the file, in bytes
    uint64_t header_size_ = 0;
    /// size of a single frame on the file, in bytes
    uint64_t frame_size_ = 0;
    /// size of the first frame on the file, in bytes (this differs from
    /// `frame_size_` only for files with fixed atoms)
    uint64_t first_frame_size_ = 0;

    /// Total number of atoms in this file
    size_t n_atoms_ = 0;
    /// Total number of free (non-fixed) atoms in this file
    size_t n_free_atoms_ = 0;

    /// Additional metadata when the file contains fixed atoms
    struct fixed_atom_data_t {
        /// is this atom fixed?
        bool fixed;
        union {
            /// valid if the atom is not fixed: index where to look for the
            /// coordinate of this atom after the first frame
            size_t free_index;
            /// valid if the atom is fixed: coordinate of the atom
            Vector3D fixed_coord;
        };
    };
    /// Metadata when there are fixed atoms in the file. This will tell for each
    /// atom if it is fixed or mobile.
    ///
    /// This is empty if all atoms are mobile (which is the common case)
    std::vector<fixed_atom_data_t> fixed_atoms_;

    /// total number of frames in the file
    size_t n_frames_ = 0;
    /// simulation timestep metadata
    struct {
        double dt = 0;
        size_t start = 0;
        size_t step = 1;
    } timesteps_;
    /// title of the file
    std::string title_;

    /// next step to read
    size_t step_ = 0;

    /// temporary buffer used when reading/writing coordinates
    std::vector<float> buffer_;
};

template<> const FormatMetadata& format_metadata<DCDFormat>();

} // namespace chemfiles

#endif
