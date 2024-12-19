// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_FORMAT_AMBER_NETCDF_HPP
#define CHEMFILES_FORMAT_AMBER_NETCDF_HPP

#include <cstddef>
#include <string>
#include <vector>

#include "chemfiles/File.hpp"
#include "chemfiles/Format.hpp"
#include "chemfiles/external/optional.hpp"

#include "chemfiles/files/Netcdf3File.hpp"

namespace chemfiles {

class Frame;
class UnitCell;
class Vector3D;
class FormatMetadata;

template <class T> class span;

/// Amber NetCDF file format implementation.
class AmberNetCDFBase: public Format {
public:
    AmberNetCDFBase(std::string convention, std::string path, File::Mode mode, File::Compression compression);

    void read(Frame& frame) final;
    void read_step(size_t step, Frame& frame) final;
    void write(const Frame& frame) override;

protected:
    struct variable_scale_t {
        netcdf3::Variable* var;
        double scale;
    };

    /// read the unit cell at the current step
    UnitCell read_cell();
    /// read the values from the variable at the current internal step to the array
    void read_array(variable_scale_t& variable, span<Vector3D> array);

    /// write the unit cell at the current step
    void write_cell(const UnitCell& cell);
    /// write the values from the array to the variable at the current internal step
    void write_array(variable_scale_t& variable, span<const Vector3D> array);

    /// associated NetCDF file.
    netcdf3::Netcdf3File file_;
    /// convention used
    std::string convention_;
    /// last step read
    size_t step_;

    struct {
        variable_scale_t coordinates;
        variable_scale_t velocities;
        variable_scale_t cell_lengths;
        variable_scale_t cell_angles;
    } variables_;

    optional<std::string> file_title_;
    size_t n_atoms_;

    std::vector<float> buffer_f32_;
    std::vector<double> buffer_f64_;

    virtual void initialize(const Frame& frame) = 0;

private:
    /// Validate the common bits between AMBER and AMBERRESTART conventions
    void validate_common();

    variable_scale_t get_variable(const std::string& name);
};

/// Amber NetCDF trajectory file format
class AmberTrajectory final: public AmberNetCDFBase {
public:
    AmberTrajectory(std::string path, File::Mode mode, File::Compression compression);

    size_t nsteps() override;
    void initialize(const Frame& frame) override;

private:
    void validate();
};

/// Amber NetCDF restart file format
class AmberRestart final: public AmberNetCDFBase {
public:
    AmberRestart(std::string path, File::Mode mode, File::Compression compression);

    void write(const Frame& frame) override;
    size_t nsteps() override;
    void initialize(const Frame& frame) override;

private:
    void validate();
};

template<> const FormatMetadata& format_metadata<AmberTrajectory>();
template<> const FormatMetadata& format_metadata<AmberRestart>();

} // namespace chemfiles

#endif
