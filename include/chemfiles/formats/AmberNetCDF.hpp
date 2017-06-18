// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_FORMAT_NC_HPP
#define CHEMFILES_FORMAT_NC_HPP

#include "chemfiles/types.hpp"
#include "chemfiles/Format.hpp"
#include "chemfiles/files/NcFile.hpp"
#include "chemfiles/FormatFactory.hpp"

namespace chemfiles {

class UnitCell;
class Topology;

/// [Amber NetCDF][NetCDF] file format reader.
///
/// [NetCDF]: http://ambermd.org/netcdf/nctraj.xhtml
class AmberNetCDFFormat final: public Format {
public:
    AmberNetCDFFormat(const std::string& path, File::Mode mode);

    void read_step(size_t step, Frame& frame) override;
    void read(Frame& frame) override;
    void write(const Frame& frame) override;

    size_t nsteps() override;
    std::string description() const override;

    // Register the Amber NetCDF format with the ".nc" extension and the
    // "AmberNetCDF" description.
    FORMAT_NAME(AmberNetCDF)
    FORMAT_EXTENSION(.nc)
private:
    /// Read the unit cell at the current internal step, the file is assumed to be valid.
    UnitCell read_cell();
    /// Generic function to read an Array3D at the current internal step,
    /// the file is assumed to be valid.
    void read_array3D(Span3D array, const std::string& name);

    /// Write an Array3D to the file, as a variable with the name `name`, at
    /// the current internal step.
    void write_array3D(const Array3D& array, const std::string& name);
    /// Write an UnitCell to the file, at the current internal step
    void write_cell(const UnitCell& cell);

    /// Associated NetCDF file.
    NcFile file_;
    /// Last read step
    size_t step_;
    /// Was the associated file validated?
    bool validated_;
};

} // namespace chemfiles

#endif
