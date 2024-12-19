// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_XDR_FILE_HPP
#define CHEMFILES_XDR_FILE_HPP

#include <cstddef>
#include <cstdint>

#include <string>
#include <vector>

#include "chemfiles/File.hpp"
#include "chemfiles/UnitCell.hpp"

#include "chemfiles/files/BinaryFile.hpp"

namespace chemfiles {

/// Partial implementation of XDR according to RFC 4506
/// (see: https://datatracker.ietf.org/doc/html/rfc4506)
/// Including additional helper routines for GROMACS
class XDRFile final : public BigEndianFile {
  public:
    XDRFile(std::string path, File::Mode mode);
    XDRFile(XDRFile&& other) noexcept = default;
    XDRFile& operator=(XDRFile&& other) = default;
    ~XDRFile() noexcept override = default;

    using BinaryFile::read_f32;
    using BinaryFile::read_f64;
    using BinaryFile::read_i32;
    using BinaryFile::write_f32;

    /// Read a non-compliant GROMACS string
    /// A GROMACS string stores the length of the string including the NULL
    /// terminator as int32 before the XDR compliant string data without the
    /// terminator. An XDR string has the same representation as opaque data.
    std::string read_gmx_string();
    /// Write a non-compliant GROMACS string
    void write_gmx_string(const std::string& value);

    /// Read compressed GROMACS floats and returns the precision
    float read_gmx_compressed_floats(std::vector<float>& data);
    /// Write compressed GROMACS floats with a given precision
    void write_gmx_compressed_floats(const std::vector<float>& data, float precision);

    /// Read the GROMACS simulation box in nano meters
    UnitCell read_gmx_box(bool use_double = false);

    /// Read an unsigned size value as i32 by performing a checked conversion
    size_t read_single_size_as_i32();

  private:
    /// Read XDR variable-length opaque data
    void read_opaque(std::vector<char>& data);
    /// Write XDR variable-length opaque data
    void write_opaque(const char* data, uint32_t count);

    /// Cache allocation for compressed data (XTC)
    std::vector<char> compressed_data_;
    /// Cache allocation for intermediate buffer (XTC)
    std::vector<int32_t> intbuf_;
};

} // namespace chemfiles

#endif
