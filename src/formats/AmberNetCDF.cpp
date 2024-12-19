// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cassert>
#include <cstddef>

#include <array>
#include <string>
#include <utility>
#include <vector>

#include "chemfiles/Property.hpp"
#include "chemfiles/config.h"
#include "chemfiles/types.hpp"
#include "chemfiles/unreachable.hpp"
#include "chemfiles/utils.hpp"
#include "chemfiles/warnings.hpp"
#include "chemfiles/error_fmt.hpp"
#include "chemfiles/external/span.hpp"
#include "chemfiles/external/optional.hpp"


#include "chemfiles/File.hpp"
#include "chemfiles/Frame.hpp"
#include "chemfiles/UnitCell.hpp"
#include "chemfiles/FormatMetadata.hpp"

#include "chemfiles/files/Netcdf3File.hpp"
#include "chemfiles/formats/AmberNetCDF.hpp"

using namespace chemfiles;

/// Create the basic builder common to AMBER and AMBERRESTART conventions
static netcdf3::Netcdf3Builder base_builder(std::string convention, const std::string& title, size_t n_atoms);

/// Find the dimension with the given name in the file
static optional<netcdf3::Dimension&> find_dimension(netcdf3::Netcdf3File& file, const std::string& name);
/// Find the id of the dimension with the given name in the builder
static size_t get_dimension_id(const netcdf3::Netcdf3Builder& builder, const std::string& name);

static double scale_for_distance(std::string unit);
static double scale_for_velocity(std::string unit);

/******************************************************************************/

AmberNetCDFBase::AmberNetCDFBase(std::string convention, std::string path, File::Mode mode, File::Compression compression):
    file_(std::move(path), mode),
    convention_(std::move(convention)),
    step_(0)
{
    if (compression != File::DEFAULT) {
        throw format_error("compression is not supported with NetCDF format");
    }

    if (!file_.initialized()) {
        // the file will be initialized on the first write
        return;
    }

    try {
        this->validate_common();
    } catch (const Error& e) {
        throw format_error("invalid Amber NetCDF file at '{}': {}", file_.path(), e.what());
    }

    const auto& attributes = file_.attributes();
    auto title = attributes.find("title");
    if (title != attributes.end() && title->second.kind() == netcdf3::Value::STRING) {
        file_title_ = title->second.as_string();
    }

    n_atoms_ = static_cast<size_t>(find_dimension(file_, "atom")->size);

    // get the variables actually defined in the file
    variables_.coordinates = this->get_variable("coordinates");
    if (variables_.coordinates.var != nullptr) {
        auto units_attr = variables_.coordinates.var->attribute("units");
        if (units_attr && units_attr->kind() == netcdf3::Value::STRING) {
            variables_.coordinates.scale *= scale_for_distance(units_attr->as_string());
        }
    } else {
        warning("Amber NetCDF reader", "the coordinates variable is not defined in this file");
    }

    variables_.velocities = this->get_variable("velocities");
    if (variables_.velocities.var != nullptr) {
        auto units_attr = variables_.velocities.var->attribute("units");
        if (units_attr && units_attr->kind() == netcdf3::Value::STRING) {
            variables_.velocities.scale *= scale_for_velocity(units_attr->as_string());
        }
    }

    variables_.cell_lengths = this->get_variable("cell_lengths");
    variables_.cell_angles = this->get_variable("cell_angles");
    if ((variables_.cell_lengths.var != nullptr) && (variables_.cell_angles.var != nullptr)) {
        auto units_attr = variables_.cell_lengths.var->attribute("units");
        if (units_attr && units_attr->kind() == netcdf3::Value::STRING) {
            variables_.cell_lengths.scale *= scale_for_distance(units_attr->as_string());
        }

        units_attr = variables_.cell_angles.var->attribute("units");
        if (units_attr && units_attr->kind() == netcdf3::Value::STRING) {
            auto units = units_attr->as_string();
            to_ascii_lowercase(units);

            if (units == "" || units == "degrees" || units == "degree") {
                // nothing to do
            } else if (units == "radians" || units == "radians") {
                constexpr double pi = 3.141592653589793238463;
                variables_.cell_angles.scale *= 180 / pi;
            } else {
                warning("Amber NetCDF reader", "unknown unit ({}) for angles", units);
            }
        }
    } else if (variables_.cell_lengths.var != nullptr) {
        if (variables_.cell_angles.var == nullptr) {
            throw format_error(
                "invalid Amber NetCDF file at '{}': "
                "cell_lengths requires cell_angles to be defined",
                file_.path()
            );
        }

        if (variables_.coordinates.var == nullptr) {
            throw format_error(
                "invalid Amber NetCDF file at '{}': "
                "cell_lengths requires coordinates to be defined",
                file_.path()
            );
        }
    }

    if (mode == File::APPEND) {
        // start writing at the end of pre-existing files in append mode
        step_ = static_cast<size_t>(file_.n_records());
    }
}

void AmberNetCDFBase::validate_common() {
    auto conventions = file_.attribute("Conventions");
    if (!conventions) {
        throw format_error("expected a 'Conventions' attribute to be defined");
    } else if (conventions->as_string() != convention_) {
        throw format_error("expected '{}' convention", convention_);
    }

    auto convention_version = file_.attribute("ConventionVersion");
    if (!convention_version) {
        throw format_error("expected a 'ConventionVersion' attribute to be defined");
    } else if (convention_version->as_string() != "1.0") {
        throw format_error("expected version 1.0 of '{}' convention", convention_);
    }

    auto spatial = find_dimension(file_, "spatial");
    if (!spatial) {
        throw format_error("missing 'spatial' dimension");
    } else if (spatial->size != 3) {
        throw format_error("'spatial' dimension must have a size of 3, got {}", spatial->size);
    }

    if (!find_dimension(file_, "atom")) {
        throw format_error("missing 'atom' dimension");
    }

    auto cell_spatial = find_dimension(file_, "cell_spatial");
    if (cell_spatial && cell_spatial->size != 3) {
        throw format_error("'cell_spatial' dimension must have a size of 3, got {}", cell_spatial->size);
    }

    auto cell_angular = find_dimension(file_, "cell_angular");
    if (cell_angular && cell_angular->size != 3) {
        throw format_error("'cell_angular' dimension must have a size of 3, got {}", cell_angular->size);
    }
}

AmberNetCDFBase::variable_scale_t AmberNetCDFBase::get_variable(const std::string& name) {
    double scale = 1.0;
    auto variable = file_.variable(name);
    if (variable) {
        auto scale_factor = variable->attribute("scale_factor");
        if (scale_factor) {
            if (scale_factor->kind() == netcdf3::Value::FLOAT) {
                scale = static_cast<double>(scale_factor->as_f32());
            } else if (scale_factor->kind() == netcdf3::Value::DOUBLE) {
                // this is non standard, but LAMMPS does it
                scale = scale_factor->as_f64();
            } else {
                throw format_error("scale_factor attribute for '{}' must be a floating point value", name);
            }
        }

        return variable_scale_t{&variable.value(), scale};
    } else {
        return variable_scale_t{nullptr, scale};
    }
}

/******************************************************************************/

void AmberNetCDFBase::read(Frame& frame) {
    this->read_step(step_, frame);
    step_++;
}

void AmberNetCDFBase::read_step(const size_t step, Frame& frame) {
    // Set the internal step_ before further reading
    step_ = step;

    frame.set_cell(read_cell());

    if (file_title_) {
        frame.set("name", file_title_.value());
    }

    frame.resize(n_atoms_);

    if (variables_.coordinates.var != nullptr) {
        this->read_array(variables_.coordinates, frame.positions());
    }

    if (variables_.velocities.var != nullptr) {
        frame.add_velocities();
        this->read_array(variables_.velocities, *frame.velocities());
    }
}

void AmberNetCDFBase::write(const Frame& frame) {
    if (!file_.initialized()) {
        this->initialize(frame);

        file_.variable("spatial")->write(0, "xyz", 3);
        file_.variable("cell_spatial")->write(0, "abc", 3);
        file_.variable("cell_angular")->write(0, "alphabeta gamma", 15);

        variables_.coordinates = this->get_variable("coordinates");
        variables_.velocities = this->get_variable("velocities");
        variables_.cell_lengths = this->get_variable("cell_lengths");
        variables_.cell_angles = this->get_variable("cell_angles");

        n_atoms_ = frame.size();
    }

    file_.add_record();

    if (frame.size() != n_atoms_) {
        throw format_error(
            "this file can only write frame with {} atoms, but the frame contains {} atoms",
            n_atoms_, frame.size()
        );
    }

    write_cell(frame.cell());

    if (variables_.coordinates.var != nullptr) {
        this->write_array(variables_.coordinates, frame.positions());
    }

    if (frame.velocities()) {
        if (variables_.velocities.var != nullptr) {
            this->write_array(variables_.velocities, *frame.velocities());
        } else {
            warning("AMBER NetCDF", "this file does not contain space for velocities, they will not be saved");
        }
    }

    step_ += 1;
}

/******************************************************************************/

UnitCell AmberNetCDFBase::read_cell() {
    if ((variables_.cell_lengths.var == nullptr) || (variables_.cell_angles.var == nullptr)) {
        // No cell information
        return {};
    }

    auto data_f32 = std::array<float, 3>();
    auto data_f64 = std::array<double, 3>();

    Vector3D lengths;
    auto& cell_lengths = variables_.cell_lengths.var;
    if (cell_lengths->type() == netcdf3::constants::NC_FLOAT) {
        cell_lengths->read(step_, data_f32.data(), data_f32.size());
        lengths = Vector3D(
            static_cast<double>(data_f32[0]),
            static_cast<double>(data_f32[1]),
            static_cast<double>(data_f32[2])
        );
    } else if (cell_lengths->type() == netcdf3::constants::NC_DOUBLE) {
        cell_lengths->read(step_, data_f64.data(), data_f64.size());
        lengths = Vector3D(
            data_f64[0],
            data_f64[1],
            data_f64[2]
        );
    } else {
        unreachable();
    }

    Vector3D angles;
    auto& cell_angles = variables_.cell_angles.var;
    if (cell_angles->type() == netcdf3::constants::NC_FLOAT) {
        cell_angles->read(step_, data_f32.data(), data_f32.size());
        angles = Vector3D(
            static_cast<double>(data_f32[0]),
            static_cast<double>(data_f32[1]),
            static_cast<double>(data_f32[2])
        );
    } else if (cell_angles->type() == netcdf3::constants::NC_DOUBLE) {
        cell_angles->read(step_, data_f64.data(), data_f64.size());
        angles = Vector3D(
            data_f64[0],
            data_f64[1],
            data_f64[2]
        );
    } else {
        unreachable();
    }

    return UnitCell(
        variables_.cell_lengths.scale * lengths,
        variables_.cell_angles.scale * angles
    );
}

void AmberNetCDFBase::read_array(variable_scale_t& variable, span<Vector3D> array) {
    if (variable.var->type() == netcdf3::constants::NC_FLOAT) {
        variable.var->read(step_, buffer_f32_);
        for (size_t i=0; i<n_atoms_; i++) {
            array[i][0] = variable.scale * static_cast<double>(buffer_f32_[3 * i + 0]);
            array[i][1] = variable.scale * static_cast<double>(buffer_f32_[3 * i + 1]);
            array[i][2] = variable.scale * static_cast<double>(buffer_f32_[3 * i + 2]);
        }
    } else if (variable.var->type() == netcdf3::constants::NC_DOUBLE) {
        variable.var->read(step_, buffer_f64_);
        for (size_t i=0; i<n_atoms_; i++) {
            array[i][0] = variable.scale * buffer_f64_[3 * i + 0];
            array[i][1] = variable.scale * buffer_f64_[3 * i + 1];
            array[i][2] = variable.scale * buffer_f64_[3 * i + 2];
        }
    } else {
        throw format_error("invalid type for variable, expected floating point");
    }
}

/******************************************************************************/

void AmberNetCDFBase::write_cell(const UnitCell& cell) {
    if ((variables_.cell_lengths.var == nullptr) || (variables_.cell_angles.var == nullptr)) {
        // no cell information
        if (cell.shape() != UnitCell::INFINITE) {
            warning("AMBER NetCDF", "this file does not contain space for unit cell data, it will not be saved");
        }
        return;
    }

    auto lengths = cell.lengths();
    auto& cell_lengths = variables_.cell_lengths.var;
    if (cell_lengths->type() == netcdf3::constants::NC_FLOAT) {
        auto data_f32 = std::array<float, 3>{
            static_cast<float>(lengths[0]),
            static_cast<float>(lengths[1]),
            static_cast<float>(lengths[2]),
        };
        cell_lengths->write(step_, data_f32.data(), data_f32.size());
    } else if (cell_lengths->type() == netcdf3::constants::NC_DOUBLE) {
        cell_lengths->write(step_, lengths.data(), 3);
    } else {
        unreachable();
    }

    auto angles = cell.angles();
    auto& cell_angles = variables_.cell_angles.var;
    if (cell_angles->type() == netcdf3::constants::NC_FLOAT) {
        auto data_f32 = std::array<float, 3>{
            static_cast<float>(angles[0]),
            static_cast<float>(angles[1]),
            static_cast<float>(angles[2]),
        };
        cell_angles->write(step_, data_f32.data(), data_f32.size());
    } else if (cell_angles->type() == netcdf3::constants::NC_DOUBLE) {
        cell_angles->write(step_, angles.data(), 3);
    } else {
        unreachable();
    }
}

void AmberNetCDFBase::write_array(variable_scale_t& variable, span<const Vector3D> array) {
    if (variable.var->type() == netcdf3::constants::NC_FLOAT) {
        buffer_f32_.resize(3 * array.size());
        for (size_t i=0; i<n_atoms_; i++) {
            buffer_f32_[3 * i + 0] = static_cast<float>(array[i][0]);
            buffer_f32_[3 * i + 1] = static_cast<float>(array[i][1]);
            buffer_f32_[3 * i + 2] = static_cast<float>(array[i][2]);
        }
        variable.var->write(step_, buffer_f32_);
    } else if (variable.var->type() == netcdf3::constants::NC_DOUBLE) {
        variable.var->write(step_, array[0].data(), 3 * array.size());
    } else {
        throw format_error("invalid type for variable, expected floating point");
    }
}

/******************************************************************************/

AmberTrajectory::AmberTrajectory(std::string path, File::Mode mode, File::Compression compression):
    AmberNetCDFBase("AMBER", std::move(path), mode, compression)
{
    if (!file_.initialized()) {
        // skip validation, the file will be initialized later
        return;
    }

    try {
        this->validate();
    } catch (const Error& e) {
        throw format_error("invalid Amber NetCDF trajectory at '{}': {}", file_.path(), e.what());
    }
}

size_t AmberTrajectory::nsteps() {
    return static_cast<size_t>(file_.n_records());
}

void AmberTrajectory::validate() {
    auto frame = find_dimension(file_, "frame");
    if (!frame) {
        throw format_error("missing 'frame' dimension");
    } else if (!frame->is_record()) {
        throw format_error("the 'frame' dimension must be the record dimension");
    }

    auto variables = file_.variables();

    auto validate_atom_variable = [](std::string name, const netcdf3::Variable& variable) {
        if (variable.type() != netcdf3::constants::NC_DOUBLE && variable.type() != netcdf3::constants::NC_FLOAT) {
            throw format_error(
                "'{}' variable must contain floating point data, got netcdf type {} instead",
                name, variable.type()
            );
        }

        auto dimensions = variable.dimensions();
        if (dimensions.size() != 3) {
            throw format_error("'{}' variable must have three dimensions", name);
        } else if (dimensions[0]->name != "frame") {
            throw format_error("first dimension of '{}' variable must be 'frame'", name);
        } else if (dimensions[1]->name != "atom") {
            throw format_error("second dimension of '{}' variable must be 'atom'", name);
        } else if (dimensions[2]->name != "spatial") {
            throw format_error("third dimension of '{}' variable must be 'spatial'", name);
        }
    };

    auto coordinates = variables.find("coordinates");
    if (coordinates != variables.end()) {
        validate_atom_variable("coordinates", coordinates->second);
    }

    auto velocities = variables.find("velocities");
    if (velocities != variables.end()) {
        validate_atom_variable("velocities", velocities->second);
    }

    auto validate_cell_variable = [](std::string name, std::string dimension, const netcdf3::Variable& variable) {
        if (variable.type() != netcdf3::constants::NC_DOUBLE && variable.type() != netcdf3::constants::NC_FLOAT) {
            throw format_error(
                "'{}' variable must contain floating point data, got netcdf type {} instead",
                name, variable.type()
            );
        }

        auto dimensions = variable.dimensions();
        if (dimensions.size() != 2) {
            throw format_error("'{}' variable must have two dimensions", name);
        } else if (dimensions[0]->name != "frame") {
            throw format_error("first dimension of '{}' variable must be 'frame'", name);
        } else if (dimensions[1]->name != dimension) {
            throw format_error("second dimension of '{}' variable must be '{}'", name, dimension);
        }
    };

    auto cell_lengths = variables.find("cell_lengths");
    if (cell_lengths != variables.end()) {
        validate_cell_variable("cell_lengths", "cell_spatial", cell_lengths->second);
    }

    auto cell_angles = variables.find("cell_angles");
    if (cell_angles != variables.end()) {
        validate_cell_variable("cell_angles", "cell_angular", cell_angles->second);
    }
}

void AmberTrajectory::initialize(const Frame& frame) {
    netcdf3::Netcdf3Builder builder = base_builder(
        "AMBER",
        frame.get<Property::STRING>("name").value_or(""),
        frame.size()
    );

    auto frame_dim = builder.add_dimension("frame", 0);
    auto atom_dim = get_dimension_id(builder, "atom");
    auto spatial_dim = get_dimension_id(builder, "spatial");
    auto cell_spatial_dim = get_dimension_id(builder, "cell_spatial");
    auto cell_angular_dim = get_dimension_id(builder, "cell_angular");

    builder.add_variable("coordinates", {
        /* type = */ netcdf3::constants::NC_FLOAT,
        /* dimensions = */ {frame_dim, atom_dim, spatial_dim},
        /* attributes = */ {{"units", netcdf3::Value("angstrom")}}
    });

    builder.add_variable("cell_lengths", {
        /* type = */ netcdf3::constants::NC_FLOAT,
        /* dimensions = */ {frame_dim, cell_spatial_dim},
        /* attributes = */ {{"units", netcdf3::Value("angstrom")}}
    });

    builder.add_variable("cell_angles", {
        /* type = */ netcdf3::constants::NC_FLOAT,
        /* dimensions = */ {frame_dim, cell_angular_dim},
        /* attributes = */ {{"units", netcdf3::Value("degree")}}
    });

    if (frame.velocities()) {
        builder.add_variable("velocities", {
            /* type = */ netcdf3::constants::NC_FLOAT,
            /* dimensions = */ {frame_dim, atom_dim, spatial_dim},
            /* attributes = */ {{"units", netcdf3::Value("angstrom/picosecond")}},
        });
    }

    std::move(builder).initialize(&file_);
}

/******************************************************************************/

AmberRestart::AmberRestart(std::string path, File::Mode mode, File::Compression compression):
    AmberNetCDFBase("AMBERRESTART", std::move(path), mode, compression)
{
    try {
        this->validate();
    } catch (const Error& e) {
        throw format_error("invalid Amber NetCDF restart at '{}': {}", file_.path(), e.what());
    }
}

void AmberRestart::write(const Frame& frame) {
    if (step_ != 0) {
        throw format_error("AMBER Restart format only supports writing one frame");
    }
    AmberNetCDFBase::write(frame);
}

size_t AmberRestart::nsteps() {
    return 1;
}

void AmberRestart::validate() {
    auto variables = file_.variables();
    auto validate_atom_variable = [](std::string name, const netcdf3::Variable& variable) {
        if (variable.type() != netcdf3::constants::NC_DOUBLE && variable.type() != netcdf3::constants::NC_FLOAT) {
            throw format_error(
                "'{}' variable must contain floating point data, got netcdf type {} instead",
                name, variable.type()
            );
        }

        auto dimensions = variable.dimensions();
        if (dimensions.size() != 2) {
            throw format_error("'{}' variable must have two dimensions", name);
        } else if (dimensions[0]->name != "atom") {
            throw format_error("first dimension of '{}' variable must be 'atom'", name);
        } else if (dimensions[1]->name != "spatial") {
            throw format_error("second dimension of '{}' variable must be 'spatial'", name);
        }
    };

    auto coordinates = variables.find("coordinates");
    if (coordinates != variables.end()) {
        validate_atom_variable("coordinates", coordinates->second);
    }

    auto velocities = variables.find("velocities");
    if (velocities != variables.end()) {
        validate_atom_variable("velocities", velocities->second);
    }

    auto validate_cell_variable = [](std::string name, std::string dimension, const netcdf3::Variable& variable) {
        if (variable.type() != netcdf3::constants::NC_DOUBLE && variable.type() != netcdf3::constants::NC_FLOAT) {
            throw format_error(
                "'{}' variable must contain floating point data, got netcdf type {} instead",
                name, variable.type()
            );
        }

        auto dimensions = variable.dimensions();
        if (dimensions.size() != 1) {
            throw format_error("'{}' variable must have one dimension", name);
        } else if (dimensions[0]->name != dimension) {
            throw format_error("first dimension of '{}' variable must be '{}'", name, dimension);
        }
    };

    auto cell_lengths = variables.find("cell_lengths");
    if (cell_lengths != variables.end()) {
        validate_cell_variable("cell_lengths", "cell_spatial", cell_lengths->second);
    }

    auto cell_angles = variables.find("cell_angles");
    if (cell_angles != variables.end()) {
        validate_cell_variable("cell_angles", "cell_angular", cell_angles->second);
    }
}

void AmberRestart::initialize(const Frame& frame) {
    netcdf3::Netcdf3Builder builder = base_builder(
        "AMBERRESTART",
        frame.get<Property::STRING>("name").value_or(""),
        frame.size()
    );

    auto atom_dim = get_dimension_id(builder, "atom");
    auto spatial_dim = get_dimension_id(builder, "spatial");
    auto cell_spatial_dim = get_dimension_id(builder, "cell_spatial");
    auto cell_angular_dim = get_dimension_id(builder, "cell_angular");

    builder.add_variable("coordinates", {
        /* type = */ netcdf3::constants::NC_DOUBLE,
        /* dimensions = */ {atom_dim, spatial_dim},
        /* attributes = */ {{"units", netcdf3::Value("angstrom")}}
    });

    builder.add_variable("cell_lengths", {
        /* type = */ netcdf3::constants::NC_DOUBLE,
        /* dimensions = */ {cell_spatial_dim},
        /* attributes = */ {{"units", netcdf3::Value("angstrom")}}
    });

    builder.add_variable("cell_angles", {
        /* type = */ netcdf3::constants::NC_DOUBLE,
        /* dimensions = */ {cell_angular_dim},
        /* attributes = */ {{"units", netcdf3::Value("degree")}}
    });

    if (frame.velocities()) {
        builder.add_variable("velocities", {
            /* type = */ netcdf3::constants::NC_DOUBLE,
            /* dimensions = */ {atom_dim, spatial_dim},
            /* attributes = */ {{"units", netcdf3::Value("angstrom/femptosecond")}}
        });
    }

    std::move(builder).initialize(&file_);
}

/******************************************************************************/

template<> const FormatMetadata& chemfiles::format_metadata<AmberTrajectory>() {
    static FormatMetadata metadata;
    metadata.name = "Amber NetCDF";
    metadata.extension = ".nc";
    metadata.description = "Amber convention for binary NetCDF molecular trajectories";
    metadata.reference = "http://ambermd.org/netcdf/nctraj.xhtml";

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

template<> const FormatMetadata& chemfiles::format_metadata<AmberRestart>() {
    static FormatMetadata metadata;
    metadata.name = "Amber Restart";
    metadata.extension = ".ncrst";
    metadata.description = "Amber convention for binary NetCDF restart files";
    metadata.reference = "http://ambermd.org/netcdf/nctraj.xhtml";

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

/******************************************************************************/

optional<netcdf3::Dimension&> find_dimension(netcdf3::Netcdf3File& file, const std::string& name) {
    for (const auto& dimension: file.dimensions()) {
        if (dimension->name == name) {
            return *dimension;
        }
    }
    return nullopt;
}

size_t get_dimension_id(const netcdf3::Netcdf3Builder& builder, const std::string& name) {
    const auto& dimensions = builder.dimensions();
    for (size_t i=0; i<dimensions.size(); i++) {
        if (dimensions[i]->name == name) {
            return i;
        }
    }
    throw error("internal error: unable to find a dimension named {} in this Netcdf3Builder", name);
}

netcdf3::Netcdf3Builder base_builder(std::string convention, const std::string& title, size_t n_atoms) {
    auto builder = netcdf3::Netcdf3Builder();
    builder.add_attribute("Conventions", {std::move(convention)});
    builder.add_attribute("ConventionVersion", {"1.0"});
    builder.add_attribute("program", {"chemfiles"});
    builder.add_attribute("programVersion", {CHEMFILES_VERSION});

    if (!title.empty()) {
        builder.add_attribute("title", {title});
    }

    auto spatial_dim = builder.add_dimension("spatial", 3);
    builder.add_dimension("atom", n_atoms);

    auto cell_spatial_dim = builder.add_dimension("cell_spatial", 3);
    auto cell_angular_dim = builder.add_dimension("cell_angular", 3);
    auto label_dim = builder.add_dimension("label", 5);

    builder.add_variable("spatial", {
        /* type = */ netcdf3::constants::NC_CHAR,
        /* dimensions = */ {spatial_dim},
        /* attributes = */ {}
    });

    builder.add_variable("cell_spatial", {
        /* type = */ netcdf3::constants::NC_CHAR,
        /* dimensions = */ {cell_spatial_dim},
        /* attributes = */ {}
    });

    builder.add_variable("cell_angular", {
        /* type = */ netcdf3::constants::NC_CHAR,
        /* dimensions = */ {cell_angular_dim, label_dim},
        /* attributes = */ {}
    });

    return builder;
}

double scale_for_distance(std::string units) {
    to_ascii_lowercase(units);

    if (units == "angstroms" || units == "angstrom" || units == "a") {
        return 1.0;
    } else if (units == "meters" || units == "meter" || units == "m") {
        return 1e10;
    } else if (units == "centimeters" || units == "centimeter" || units == "cm") {
        return 1e8;
    } else if (units == "micrometers" || units == "micrometer" || units == "µm"  || units == "um") {
        return 1e4;
    } else if (units == "nanometers" || units == "nanometer" || units == "nm") {
        return 1e4;
    } else if (units == "bohrs" || units == "bohr") {
        return 0.52918;
    } else {
        warning("Amber NetCDF reader", "unknown unit ({}) for distances", units);
        return 1.0;
    }
}

double scale_for_velocity(std::string units) {
    to_ascii_lowercase(units);

    auto splitted = split(units, '/');
    if (splitted.size() != 2) {
        warning("Amber NetCDF reader", "unknown unit ({}) for velocities", units);
        return 1.0;
    }
    auto scale = scale_for_distance(std::string(splitted[0]));
    auto time_unit = splitted[1];

    if (time_unit == "picoseconds" || time_unit == "picosecond" || time_unit == "ps") {
        // nothing to do
    } else if (time_unit == "femtoseconds" || time_unit == "femtosecond" || time_unit == "fs") {
        scale *= 1e3;
    } else if (time_unit == "nanoseconds" || time_unit == "nanosecond" || time_unit == "ns") {
        scale *= 1e-3;
    } else if (time_unit == "microseconds" || time_unit == "microsecond" || time_unit == "µs" || time_unit == "us") {
        scale *= 1e-6;
    } else if (time_unit == "seconds" || time_unit == "second" || time_unit == "s") {
        scale *= 1e-12;
    } else {
        warning("Amber NetCDF reader", "unknown unit ({}) for time", time_unit);
    }

    return scale;
}
