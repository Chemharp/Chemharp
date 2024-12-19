// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cassert>

#include <map>
#include <array>
#include <string>
#include <vector>
#include <string_view>

#include "chemfiles/types.hpp"
#include "chemfiles/utils.hpp"
#include "chemfiles/parse.hpp"
#include "chemfiles/warnings.hpp"
#include "chemfiles/error_fmt.hpp"
#include "chemfiles/external/optional.hpp"

#include "chemfiles/File.hpp"
#include "chemfiles/Atom.hpp"
#include "chemfiles/Frame.hpp"
#include "chemfiles/Property.hpp"
#include "chemfiles/Residue.hpp"
#include "chemfiles/Topology.hpp"
#include "chemfiles/UnitCell.hpp"
#include "chemfiles/FormatMetadata.hpp"

#include "chemfiles/formats/GRO.hpp"

using namespace chemfiles;

template<> const FormatMetadata& chemfiles::format_metadata<GROFormat>() {
    static FormatMetadata metadata;
    metadata.name = "GRO";
    metadata.extension = ".gro";
    metadata.description = "GROMACS GRO text format";
    metadata.reference = "http://manual.gromacs.org/current/reference-manual/file-formats.html#gro";

    metadata.read = true;
    metadata.write = true;
    metadata.memory = true;

    metadata.positions = true;
    metadata.velocities = true;
    metadata.unit_cell = true;
    metadata.atoms = true;
    metadata.bonds = false;
    metadata.residues = true;
    return metadata;
}

using chemfiles::private_details::is_upper_triangular;

/// Check the number of digits before the decimal separator to be sure than
/// we can represent them. In case of error, use the given `context` in the error
/// message
static void check_values_size(const Vector3D& values, unsigned width, const std::string& context);

void GROFormat::read_next(Frame& frame) {
    residues_.clear();

    // GRO comment line is used as frame name
    auto frame_name = trim(file_.readline());
    if (!frame_name.empty()) {
        frame.set("name", std::string(frame_name));
    }

    size_t natoms = 0;
    try {
        natoms = parse<size_t>(file_.readline());
    } catch (const Error& e) {
        throw format_error("can not read number of atoms in GRO file: {}", e.what());
    }

    frame.add_velocities();
    frame.reserve(natoms);

    for (size_t i=0; i<natoms; i++) {
        auto line = file_.readline();
        if (line.length() < 44) {
            throw format_error("GRO Atom line is too small: '{}'", line);
        }

        optional<int64_t> resid = nullopt;
        try {
            resid = parse<int64_t>(line.substr(0, 5));
        } catch (const Error&) {
            // Invalid residue, we'll skip it
            warning("GRO Reader", "skiping invalid residue with resid '{}'", line.substr(0, 5));
        }

        auto resname = std::string(trim(line.substr(5, 5)));
        auto name = std::string(trim(line.substr(10, 5)));

        // GRO files store atoms in nanometer, we need to convert to Angstroms
        auto x = parse<double>(line.substr(20, 8)) * 10;
        auto y = parse<double>(line.substr(28, 8)) * 10;
        auto z = parse<double>(line.substr(36, 8)) * 10;

        double vx = 0;
        double vy = 0;
        double vz = 0;
        if (line.length() >= 68) {
            vx = parse<double>(line.substr(44, 8)) * 10;
            vy = parse<double>(line.substr(52, 8)) * 10;
            vz = parse<double>(line.substr(60, 8)) * 10;
        }
        frame.add_atom(Atom(name), {x, y, z}, {vx, vy, vz});

        if (!resid) {
            continue;
        }

        if (residues_.find(*resid) == residues_.end()) {
            Residue residue(resname, *resid);
            residue.add_atom(frame.size() - 1);

            residues_.insert({*resid, residue});
        } else {
            // Just add this atom to the residue
            residues_.at(*resid).add_atom(frame.size() - 1);
        }
    }

    auto box = file_.readline();
    auto box_values = split(box, ' ');

    if (box_values.size() == 3) {
        auto lengths = Vector3D(
            parse<double>(box_values[0]) * 10,
            parse<double>(box_values[1]) * 10,
            parse<double>(box_values[2]) * 10
        );

        frame.set_cell({lengths});
    } else if (box_values.size() == 9) {
        auto v1_x = parse<double>(box_values[0]) * 10;
        auto v2_y = parse<double>(box_values[1]) * 10;
        auto v3_z = parse<double>(box_values[2]) * 10;

        assert(parse<double>(box_values[3]) == 0);
        assert(parse<double>(box_values[4]) == 0);

        auto v2_x = parse<double>(box_values[5]) * 10;

        assert(parse<double>(box_values[6]) == 0);

        auto v3_x = parse<double>(box_values[7]) * 10;
        auto v3_y = parse<double>(box_values[8]) * 10;

        auto cell = UnitCell({
            v1_x, v2_x, v3_x,
            0.00, v2_y, v3_y,
            0.00, 0.00, v3_z
        });
        frame.set_cell(cell);
    }

    for (auto& residue: residues_) {
        frame.add_residue(residue.second);
    }
}

static std::string to_gro_index(uint64_t i) {
    if (i >= 99999) {
        if (i == 99999) {
            // Only warn once for this
            warning("GRO writer", "too many atoms, removing atomic id bigger than 100000");
        }
        return "*****";
    } else {
        return std::to_string(i + 1);
    }
}

void GROFormat::write_next(const Frame& frame) {
    file_.print("{}\n", frame.get<Property::STRING>("name").value_or("GRO File produced by chemfiles"));
    file_.print("{: >5d}\n", frame.size());

    // Only use numbers bigger than the biggest residue id as "resSeq" for
    // atoms without associated residue, and start generated residue id at
    // 1
    int64_t max_resid = 1;
    for (const auto& residue: frame.topology().residues()) {
        auto resid = residue.id();
        if (resid && resid.value() > max_resid) {
            max_resid = resid.value() + 1;
        }
    }

    const auto& positions = frame.positions();
    for (size_t i = 0; i < frame.size(); i++) {
        std::string resname = "XXXXX";
        std::string resid = "-1";
        auto residue = frame.topology().residue_for_atom(i);
        if (residue) {
            resname = residue->name();
            if (resname.length() > 5) {
                warning("GRO writer",
                    "residue '{}' name is too long, it will be truncated",
                    resname
                );
                resname = resname.substr(0, 5);
            }
        }

        if (residue && residue->id()) {
            auto value = residue->id().value();
            if (value <= 0) {
                warning("GRO writer", "the residue id '{}' should not be negative or zero, treating it as blank", value);
                value = max_resid++;
                if (value <= 99999) {
                    resid = std::to_string(value);
                }
            } else if (value <= 99999) {
                resid = std::to_string(value);
            } else {
                warning("GRO writer", "too many residues, removing residue id");
            }
        } else {
            // We need to manually assign a residue ID
            auto value = max_resid++;
            if (value <= 99999) {
                resid = std::to_string(value);
            }
        }

        assert(resname.length() <= 5);
        auto pos = positions[i] / 10;
        check_values_size(pos, 8, "atomic position");

        if (frame.velocities()) {
            auto vel = (*frame.velocities())[i] / 10;
            check_values_size(vel, 8, "atomic velocity");
            file_.print(
                "{: >5}{: <5}{: >5}{: >5}{:8.3f}{:8.3f}{:8.3f}{:8.4f}{:8.4f}{:8.4f}\n",
                resid, resname, frame[i].name(), to_gro_index(i), pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]
            );
        } else {
            file_.print(
                "{: >5}{: <5}{: >5}{: >5}{:8.3f}{:8.3f}{:8.3f}\n",
                resid, resname, frame[i].name(), to_gro_index(i), pos[0], pos[1], pos[2]
            );
        }
    }

    const auto& cell = frame.cell();
    // While this line is free form, we should try to print it in a pretty way that most gro parsers expect
    // This means we cannot support incredibly large cell sizes, but these are likely not practical anyway
    if (cell.shape() == UnitCell::ORTHORHOMBIC || cell.shape() == UnitCell::INFINITE) {
        auto lengths = cell.lengths() / 10;
        check_values_size(lengths, 8, "unit cell");
        // print zeros if the cell is infinite, this line is still required
        file_.print("   {:8.5f} {:8.5f} {:8.5f}\n", lengths[0], lengths[1], lengths[2]);
    } else { // Triclinic
        const auto& matrix = cell.matrix() / 10;
        if (!is_upper_triangular(matrix)) {
            throw format_error("unsupported triclinic but non upper-triangular cell matrix in GRO writer");
        }
        check_values_size(Vector3D(matrix[0][0], matrix[1][1], matrix[2][2]), 8, "unit cell");
        check_values_size(Vector3D(matrix[0][1], matrix[0][2], matrix[1][2]), 8, "unit cell");
        file_.print(
            "   {:8.5f} {:8.5f} {:8.5f} 0.0 0.0 {:8.5f} 0.0 {:8.5f} {:8.5f}\n",
            matrix[0][0], matrix[1][1], matrix[2][2], matrix[0][1], matrix[0][2], matrix[1][2]
        );
    }
}

void check_values_size(const Vector3D& values, unsigned width, const std::string& context) {
    double max_pos = std::pow(10.0, width) - 1;
    double max_neg = -std::pow(10.0, width - 1) + 1;
    if (values[0] > max_pos || values[1] > max_pos || values[2] > max_pos ||
        values[0] < max_neg || values[1] < max_neg || values[2] < max_neg) {
        throw format_error(
            "value in {} is too big for representation in GRO format", context
        );
    }
}

optional<uint64_t> GROFormat::forward() {
    auto position = file_.tellpos();
    size_t n_atoms = 0;

    // Skip the comment line
    file_.readline();

    if (file_.eof()) {
        return nullopt;
    }

    auto line = file_.readline();
    try {
        n_atoms = parse<size_t>(line);
    } catch (const Error&) {
        throw format_error(
            "could not read the number of atoms for GRO format: the line is '{}'",
            line
        );
    }

    for (size_t i=0; i<n_atoms+1; i++) {
        if (file_.eof()) {
            throw format_error(
                "not enough lines in '{}' for GRO format", file_.path()
            );
        }
        file_.readline();
    }

    return position;
}
