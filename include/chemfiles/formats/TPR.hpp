// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_TPR_FORMAT_HPP
#define CHEMFILES_TPR_FORMAT_HPP

#include <cstddef>
#include <cstdint>

#include <string>
#include <vector>

#include "chemfiles/File.hpp"
#include "chemfiles/Format.hpp"

#include "chemfiles/files/XDRFile.hpp"

namespace chemfiles {
class Frame;
class FormatMetadata;

/// GROMACS TPR file format reader.
/// The reader closely follows the original GROMACS implementation
/// and often refers to specific files in the GROMACS repository.
/// In the following, the repository path is abbreviated as
/// <GMX> := https://gitlab.com/gromacs/gromacs/-/tree/v2022.2
class TPRFormat final : public Format {
  public:
    TPRFormat(std::string path, File::Mode mode, File::Compression compression);

    void read_step(size_t step, Frame& frame) override;
    void read(Frame& frame) override;
    size_t nsteps() override;

  private:
    // Since GROMACS 2020 (TPR version 119) the way the body is deserialized changes.
    // For `FileIOXdr` see <GMX>/src/gromacs/fileio/gmxfio_xdr.cpp
    // and <GMX>/src/gromacs/fileio/gmx_internal_xdr.cpp
    // Deviations from the XDR specification:
    //  - uses `read_gmx_string()` from XDRFile
    //  - stores unsigned char as u32
    //  - stores unsigned short as u32
    //  - stores bool as i32
    // For `InMemory` see <GMX>/src/gromacs/utility/inmemoryserializer.cpp
    // Deviations from the XDR specification:
    //  - uses a different string serializer
    //  - stores bool as u8
    enum TprBodyConvention {
        FileIOXdr,
        InMemory,
    };

    /// The header contains information about the general aspect of the system.
    /// see `TpxFileHeader` in <GMX>/api/legacy/include/gromacs/fileio/tpxio.h
    struct TprHeader {
        /// Indicates if input record is present
        bool has_input_record = false;
        /// Indicates if a box is present
        bool has_box = false;
        /// Indicates if a topology is present
        bool has_topology = false;
        /// Indicates if coordinates are present
        bool has_positions = false;
        /// Indicates if velocities are present
        bool has_velocities = false;
        /// Indicates if forces are present
        /// No longer supported, but retained so old TPR can be read
        bool has_forces = false;
        /// The total number of atoms
        size_t natoms = 0;
        /// The number of temperature coupling groups
        size_t ngroups_temperature_coupling = 0;
        /// Current value of lambda
        double lambda = 0;
        /// File version
        int file_version = 0;
        /// File generation
        int file_generation = 0;
        /// If the TPR file was written in double precision
        bool use_double = false;
        /// Version-dependent serializer used for the body
        TprBodyConvention body_convention = FileIOXdr;

        /// Size of real values in bytes
        size_t sizeof_real() const {
            return this->use_double ? sizeof(double) : sizeof(float);
        }
    };

    /// Read the file header
    /// see `do_tpxheader()` in <GMX>/src/gromacs/fileio/tpxio.cpp
    void read_header();

    /// Read box and skip temperature coupling groups
    /// see `do_tpx_state_first()` in <GMX>/src/gromacs/fileio/tpxio.cpp
    void read_box(Frame& frame);

    /// Read the topology which contains atoms, residues, and bonds.
    /// Angles, Dihedrals, and Impropers are not added to the frame.
    /// see `do_tpx_mtop()` and `do_mtop` in <GMX>/src/gromacs/fileio/tpxio.cpp
    void read_topology(Frame& frame);

    // Read positions, velocities, and skip forces
    // see `do_tpx_state_second` in <GMX>/src/gromacs/fileio/tpxio.cpp
    void read_coordinates(Frame& frame);

    // Read all symbol strings which can be referenced by index.
    // see `do_symtab` in <GMX>/src/gromacs/fileio/tpxio.cpp
    std::vector<std::string> read_symbol_table();

    /// Read an index to an entry from the symbol table
    /// Return a reference to the entry in the table
    /// see `do_symstr()` in <GMX>/src/gromacs/fileio/tpxio.cpp
    const std::string& read_symbol_table_entry(const std::vector<std::string>& table);

    /// Read a GROMACS string according to the body convention
    std::string read_gmx_string();

    /// Read a GROMACS unsigned char according to the body convention
    uint8_t read_gmx_uchar();

    /// Read a GROMACS bool according to the body convention
    bool read_gmx_bool();

    /// Associated XDR file
    XDRFile file_;
    /// TPR header of the file
    TprHeader header_;
    /// The next step to read
    size_t step_ = 0;
};

template <> const FormatMetadata& format_metadata<TPRFormat>();

} // namespace chemfiles

#endif
