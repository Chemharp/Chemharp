// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>

#include <array>
#include <string>
#include <vector>
#include <string_view>

#include "chemfiles/Property.hpp"
#include "chemfiles/error_fmt.hpp"
#include "chemfiles/external/optional.hpp"
#include "chemfiles/parse.hpp"
#include "chemfiles/types.hpp"
#include "chemfiles/unreachable.hpp"
#include "chemfiles/utils.hpp"
#include "chemfiles/warnings.hpp"

#include "chemfiles/File.hpp"
#include "chemfiles/FormatMetadata.hpp"
#include "chemfiles/Frame.hpp"
#include "chemfiles/UnitCell.hpp"

#include "chemfiles/formats/LAMMPSTrajectory.hpp"

using namespace chemfiles;

template <> const FormatMetadata& chemfiles::format_metadata<LAMMPSTrajectoryFormat>() {
    static FormatMetadata metadata;
    metadata.name = "LAMMPS";
    metadata.extension = ".lammpstrj";
    metadata.description = "LAMMPS text trajectory format";
    metadata.reference = "https://lammps.sandia.gov/doc/dump.html";

    metadata.read = true;
    metadata.write = true;
    metadata.memory = true;

    metadata.positions = true;
    metadata.velocities = true;
    metadata.unit_cell = true;
    metadata.atoms = true;
    metadata.bonds = false;
    metadata.residues = false;
    return metadata;
}

using chemfiles::private_details::is_upper_triangular;

static optional<std::string_view> get_item(std::string_view line) {
    auto splitted = split(line, ':');
    if (splitted.size() != 2 || trim(splitted[0]) != "ITEM") {
        return nullopt;
    }
    return trim(splitted[1]);
}

std::array<double, 3> LAMMPSTrajectoryFormat::read_cell(Frame& frame) {
    auto line = file_.readline();
    auto item = get_item(line);
    try {
        if (!item) {
            throw format_error("expected an ITEM entry in LAMMPS format, got '{}'", line);
        }
        auto splitted = split(*item, ' ');
        if (splitted[0] == "BOX" && splitted[1] == "BOUNDS") {
            auto matrix = Matrix3D::unit();
            std::array<double, 3> origin;
            auto shape = UnitCell::ORTHORHOMBIC;
            if (splitted.size() >= 5 && (*item).find("xy xz yz") != std::string_view::npos) {
                shape = UnitCell::TRICLINIC;
            }
            line = file_.readline();
            splitted = split(line, ' ');
            if ((shape == UnitCell::ORTHORHOMBIC && splitted.size() != 2) ||
                (shape == UnitCell::TRICLINIC && splitted.size() != 3)) {
                size_t expected_dims = (shape == UnitCell::ORTHORHOMBIC) ? 2 : 3;
                throw format_error(
                    "incomplete box dimensions in LAMMPS format, expected {} but got {}",
                    expected_dims, splitted.size());
            }
            double xlo = parse<double>(splitted[0]);
            double xhi = parse<double>(splitted[1]);
            matrix[0][0] = xhi - xlo;
            origin[0] = xlo;
            if (shape == UnitCell::TRICLINIC) {
                matrix[0][1] = parse<double>(splitted[2]);
            }

            line = file_.readline();
            splitted = split(line, ' ');
            if ((shape == UnitCell::ORTHORHOMBIC && splitted.size() != 2) ||
                (shape == UnitCell::TRICLINIC && splitted.size() != 3)) {
                size_t expected_dims = (shape == UnitCell::ORTHORHOMBIC) ? 2 : 3;
                throw format_error(
                    "incomplete box dimensions in LAMMPS format, expected {} but got {}",
                    expected_dims, splitted.size());
            }
            double ylo = parse<double>(splitted[0]);
            double yhi = parse<double>(splitted[1]);
            matrix[1][1] = yhi - ylo;
            origin[1] = ylo;
            if (shape == UnitCell::TRICLINIC) {
                matrix[0][2] = parse<double>(splitted[2]);
            }

            line = file_.readline();
            splitted = split(line, ' ');
            if ((shape == UnitCell::ORTHORHOMBIC && splitted.size() != 2) ||
                (shape == UnitCell::TRICLINIC && splitted.size() != 3)) {
                size_t expected_dims = (shape == UnitCell::ORTHORHOMBIC) ? 2 : 3;
                throw format_error(
                    "incomplete box dimensions in LAMMPS format, expected {} but got {}",
                    expected_dims, splitted.size());
            }
            double zlo = parse<double>(splitted[0]);
            double zhi = parse<double>(splitted[1]);
            matrix[2][2] = zhi - zlo;
            origin[2] = zlo;
            if (shape == UnitCell::TRICLINIC) {
                matrix[1][2] = parse<double>(splitted[2]);
            }

            auto cell = UnitCell(matrix);
            frame.set_cell(cell);
            return origin;
        } else {
            throw format_error("missing 'BOX BOUNDS' item in LAMMPS format");
        }
    } catch (const Error& e) {
        throw format_error("can not read box header in LAMMPS format: {}", e.what());
    }
}

/// LAMMPS is able to dump various per-atom properties and arbitrary user-defined
/// variables
enum lammps_atom_attr_t {
    // other possible attributes that are not important for chemfiles
    CUSTOM,
    // atom ID
    ATOMID,
    // atom type
    TYPE,
    // name of atom element
    ELEMENT,
    // atom mass
    MASS,
    // unscaled atom coordinates
    POSX,
    POSY,
    POSZ,
    // scaled atom coordinates
    POSXS,
    POSYS,
    POSZS,
    // unwrapped atom coordinates
    POSXU,
    POSYU,
    POSZU,
    // scaled unwrapped atom coordinates
    POSXSU,
    POSYSU,
    POSZSU,
    // box image that the atom is in
    IMGX,
    IMGY,
    IMGZ,
    // atom velocities
    VELX,
    VELY,
    VELZ,
    // atom charge
    CHARGE,
};

// LAMMPS is able to dump the atomic positions in multiple formats
// multiple of these representations might be present simultaneously in a frame
// possible representations of the atomic positions
enum lammps_position_representation_t {
    // no atomic positions in frame
    NOPOS,
    // wrapped positions as x,y,z
    WRAPPED,
    // scaled positions as xs,ys,zs
    SCALED,
    // unwrapped positions as xu,yu,zu
    UNWRAPPED,
    // scaled unwrapped positions as xsu,ysu,zsu
    SCALED_UNWRAPPED,
};

static lammps_atom_attr_t attribute_from_str(std::string_view attr_str) {
    if (attr_str == "id") {
        return ATOMID;
    } else if (attr_str == "type") {
        return TYPE;
    } else if (attr_str == "element") {
        return ELEMENT;
    } else if (attr_str == "mass") {
        return MASS;
    } else if (attr_str == "x") {
        return POSX;
    } else if (attr_str == "y") {
        return POSY;
    } else if (attr_str == "z") {
        return POSZ;
    } else if (attr_str == "xs") {
        return POSXS;
    } else if (attr_str == "ys") {
        return POSYS;
    } else if (attr_str == "zs") {
        return POSZS;
    } else if (attr_str == "xu") {
        return POSXU;
    } else if (attr_str == "yu") {
        return POSYU;
    } else if (attr_str == "zu") {
        return POSZU;
    } else if (attr_str == "xsu") {
        return POSXSU;
    } else if (attr_str == "ysu") {
        return POSYSU;
    } else if (attr_str == "zsu") {
        return POSZSU;
    } else if (attr_str == "ix") {
        return IMGX;
    } else if (attr_str == "iy") {
        return IMGY;
    } else if (attr_str == "iz") {
        return IMGZ;
    } else if (attr_str == "vx") {
        return VELX;
    } else if (attr_str == "vy") {
        return VELY;
    } else if (attr_str == "vz") {
        return VELZ;
    } else if (attr_str == "q") {
        return CHARGE;
    } else {
        return CUSTOM;
    }
}

struct AtomField {
    std::string name;
    lammps_atom_attr_t kind;
};

static lammps_position_representation_t
detect_best_pos_representation(const std::vector<AtomField>& fields) {
    int wrapped_count = 0;
    int scaled_count = 0;
    int unwrapped_count = 0;
    int scaled_unwrapped_count = 0;
    for (const auto& atom_field : fields) {
        switch (atom_field.kind) {
        case POSX:
        case POSY:
        case POSZ:
            ++wrapped_count;
            break;
        case POSXS:
        case POSYS:
        case POSZS:
            ++scaled_count;
            break;
        case POSXU:
        case POSYU:
        case POSZU:
            ++unwrapped_count;
            break;
        case POSXSU:
        case POSYSU:
        case POSZSU:
            ++scaled_unwrapped_count;
            break;
        default:
            break;
        }
    }
    // choose the coordinate style that is defined for the most dimensions
    // prefer unwrapped coordinates over wrapped coordinates
    // and prefer non-scaled coordinates over scaled coordinates
    if (wrapped_count == 0 && scaled_count == 0 && unwrapped_count == 0 &&
        scaled_unwrapped_count == 0) {
        // no positions defined
        return NOPOS;
    } else if (unwrapped_count >= scaled_unwrapped_count && unwrapped_count >= wrapped_count &&
               unwrapped_count >= scaled_count) {
        return UNWRAPPED;
    } else if (scaled_unwrapped_count >= wrapped_count && scaled_unwrapped_count >= scaled_count) {
        return SCALED_UNWRAPPED;
    } else if (wrapped_count >= scaled_count) {
        return WRAPPED;
    } else {
        return SCALED;
    }
}

static void unwrap(Vector3D& position, std::array<int, 3>& image, Matrix3D& matrix) {
    // unwrap coordinates by using image data
    position[0] += image[0] * matrix[0][0] + image[1] * matrix[0][1] + image[2] * matrix[0][2];
    position[1] += image[1] * matrix[1][1] + image[2] * matrix[1][2];
    position[2] += image[2] * matrix[2][2];
}

void LAMMPSTrajectoryFormat::read_next(Frame& frame) {
    auto item = get_item(file_.readline());
    if (!item) {
        throw format_error("can not read next step as LAMMPS format: expected an ITEM entry");
    }
    if (*item == "UNITS") { // optional
        frame.set("lammps_units", std::string(trim(file_.readline())));
        item = get_item(file_.readline());
        if (!item) {
            throw format_error("can not read next step as LAMMPS format: expected an ITEM entry");
        }
    }
    if (*item == "TIME") { // optional
        double time = parse<double>(trim(file_.readline()));
        frame.set("time", time);
        item = get_item(file_.readline());
        if (!item) {
            throw format_error("can not read next step as LAMMPS format: expected an ITEM entry");
        }
    }

    if (*item == "TIMESTEP") {
        int64_t timestep = parse<int64_t>(trim(file_.readline()));
        frame.set_step(static_cast<size_t>(timestep));
    } else {
        throw format_error("can not read next step as LAMMPS format: expected 'TIMESTEP' got '{}'",
                           *item);
    }

    size_t natoms = 0;
    item = get_item(file_.readline());
    if (item && *item == "NUMBER OF ATOMS") {
        int64_t lmp_natoms = parse<int64_t>(trim(file_.readline()));
        natoms = static_cast<size_t>(lmp_natoms);
    } else {
        throw format_error("can not read next step as LAMMPS format: expected "
                           "'NUMBER OF ATOMS' got '{}'",
                           *item);
    }

    // LAMMPS can have boxes that do not use (0,0,0) as origin
    auto origin = read_cell(frame);

    item = get_item(file_.readline());
    if (!item) {
        throw format_error("can not read next step as LAMMPS format: expected an ITEM entry");
    }
    auto atoms_item = split(*item, ' ');
    if (atoms_item.empty() || atoms_item[0] != "ATOMS") {
        throw format_error("can not read next step as LAMMPS format: expected 'ATOMS' got '{}'",
                           *item);
    }
    std::vector<AtomField> fields;
    fields.reserve(atoms_item.size() - 1);
    optional<size_t> atomid_column = nullopt;
    std::vector<bool> duplicate_check;
    optional<std::vector<std::array<int, 3>>> images = nullopt;
    for (size_t i = 1; i < atoms_item.size(); ++i) {
        auto attr = attribute_from_str(atoms_item[i]);
        if (attr == ATOMID) {
            atomid_column = i - 1;
            duplicate_check = std::vector<bool>(natoms, false);
        }
        if (attr == VELX || attr == VELY || attr == VELZ) {
            frame.add_velocities();
        }
        if (attr == IMGX || attr == IMGY || attr == IMGZ) {
            images = std::vector<std::array<int, 3>>(natoms, {0, 0, 0});
        }
        fields.push_back({std::string(atoms_item[i]), attr});
    }
    lammps_position_representation_t use_pos_repr = detect_best_pos_representation(fields);

    frame.resize(natoms);
    auto positions = frame.positions();
    auto velocities = frame.velocities();

    for (size_t i = 0; i < natoms; ++i) {
        auto line = file_.readline();
        auto splitted = split(line, ' ');
        if (splitted.size() != fields.size()) {
            throw format_error(
                "LAMMPS atom line has wrong number of fields: expected {} got {}",
                fields.size(), splitted.size()
            );
        }

        size_t atomid = i;
        if (atomid_column) {
            // LAMMPS uses atom IDs that start with 1
            atomid = parse<size_t>(splitted[*atomid_column]);
            assert(atomid > 0);
            --atomid; // the frame uses zero-based indices
            assert(duplicate_check.size() == natoms);
            if (duplicate_check[atomid]) {
                throw format_error(
                    "found atoms with the same ID in LAMMPS format: {} is already present",
                    atomid + 1);
            }
            duplicate_check[atomid] = true;
        }

        auto& atom = frame[atomid];
        for (size_t j = 0; j < fields.size(); ++j) {
            switch (fields[j].kind) {
            case TYPE:
                atom.set_type(std::string(splitted[j]));
                break;
            case ELEMENT:
                atom.set_name(std::string(splitted[j]));
                break;
            case MASS:
                atom.set_mass(parse<double>(splitted[j]));
                break;
            case POSX:
                if (use_pos_repr == WRAPPED) {
                    positions[atomid][0] = parse<double>(splitted[j]);
                }
                break;
            case POSY:
                if (use_pos_repr == WRAPPED) {
                    positions[atomid][1] = parse<double>(splitted[j]);
                }
                break;
            case POSZ:
                if (use_pos_repr == WRAPPED) {
                    positions[atomid][2] = parse<double>(splitted[j]);
                }
                break;
            case POSXS:
                if (use_pos_repr == SCALED) {
                    // store scaled position (same for POSYS and POSZS)
                    // transform at the end when all three coordinates are known
                    positions[atomid][0] = parse<double>(splitted[j]);
                }
                break;
            case POSYS:
                if (use_pos_repr == SCALED) {
                    positions[atomid][1] = parse<double>(splitted[j]);
                }
                break;
            case POSZS:
                if (use_pos_repr == SCALED) {
                    positions[atomid][2] = parse<double>(splitted[j]);
                }
                break;
            case POSXU:
                if (use_pos_repr == UNWRAPPED) {
                    positions[atomid][0] = parse<double>(splitted[j]);
                }
                break;
            case POSYU:
                if (use_pos_repr == UNWRAPPED) {
                    positions[atomid][1] = parse<double>(splitted[j]);
                }
                break;
            case POSZU:
                if (use_pos_repr == UNWRAPPED) {
                    positions[atomid][2] = parse<double>(splitted[j]);
                }
                break;
            case POSXSU:
                if (use_pos_repr == SCALED_UNWRAPPED) {
                    // store scaled position (same for POSYSU and POSZSU)
                    // transform at the end when all three coordinates are known
                    positions[atomid][0] = parse<double>(splitted[j]);
                }
                break;
            case POSYSU:
                if (use_pos_repr == SCALED_UNWRAPPED) {
                    positions[atomid][1] = parse<double>(splitted[j]);
                }
                break;
            case POSZSU:
                if (use_pos_repr == SCALED_UNWRAPPED) {
                    positions[atomid][2] = parse<double>(splitted[j]);
                }
                break;
            case IMGX:
                assert(images);
                (*images)[atomid][0] = parse<int>(splitted[j]);
                break;
            case IMGY:
                assert(images);
                (*images)[atomid][1] = parse<int>(splitted[j]);
                break;
            case IMGZ:
                assert(images);
                (*images)[atomid][2] = parse<int>(splitted[j]);
                break;
            case VELX:
                assert(velocities);
                (*velocities)[atomid][0] = parse<double>(splitted[j]);
                break;
            case VELY:
                assert(velocities);
                (*velocities)[atomid][1] = parse<double>(splitted[j]);
                break;
            case VELZ:
                assert(velocities);
                (*velocities)[atomid][2] = parse<double>(splitted[j]);
                break;
            case CHARGE: {
                double charge = parse<double>(splitted[j]);
                atom.set_charge(charge);
            } break;
            case ATOMID:
                break;
            case CUSTOM:
                try {
                    // LAMMPS should always write double values
                    atom.set(fields[j].name, parse<double>(splitted[j]));
                } catch (const Error&) {
                    // use the string value as fallback
                    atom.set(fields[j].name, std::string(splitted[j]));
                }
                break;
            default:
                unreachable();
            }
        }
    }

    if (use_pos_repr == SCALED || use_pos_repr == SCALED_UNWRAPPED) {
        // all atoms currently know their scales position
        // transform the scaled coordinates to a non-scaled representation
        auto matrix = frame.cell().matrix();
        for (size_t i = 0; i < natoms; ++i) {
            // x = xlo + xs * (xhi - xlo) + ys * xy + zs * xz
            positions[i][0] = origin[0] + positions[i][0] * matrix[0][0] +
                              positions[i][1] * matrix[0][1] + positions[i][2] * matrix[0][2];
            // y = ylo + ys * (yhi - ylo) + z * yz
            positions[i][1] =
                origin[1] + positions[i][1] * matrix[1][1] + positions[i][2] * matrix[1][2];
            // z = zlo + zs * (zhi - zlo)
            positions[i][2] = origin[2] + positions[i][2] * matrix[2][2];
            if (images && use_pos_repr != SCALED_UNWRAPPED) {
                // unwrap coordinates by using image data
                unwrap(positions[i], (*images)[i], matrix);
            }
        }
    } else if (images && use_pos_repr != UNWRAPPED) {
        // unwrap coordinates by using image data
        auto matrix = frame.cell().matrix();
        for (size_t i = 0; i < natoms; ++i) {
            unwrap(positions[i], (*images)[i], matrix);
        }
    }

    if (use_pos_repr == UNWRAPPED || use_pos_repr == SCALED_UNWRAPPED || images) {
        frame.set("is_unwrapped", true);
    }
    else {
        frame.set("is_unwrapped", false);
    }
}

static optional<size_t> parse_lammps_type(const std::string& type_str) {
    if (type_str.empty()) {
        return nullopt;
}
    try {
        int type = parse<int>(type_str);
        if (type > 0) {
            return static_cast<size_t>(type);
        }
    } catch (const Error&) {
        // parsing errors indicate invalid types
    }
    return nullopt;
}

void LAMMPSTrajectoryFormat::write_next(const Frame& frame) {
    // use angstrom and femtosecond as default
    file_.print("ITEM: UNITS\n{:s}\n",
                frame.get<Property::STRING>("lammps_units").value_or("real"));
    if (frame.get("time")) {
        file_.print("ITEM: TIME\n{:.16g}\n", (*frame.get("time")).as_double());
    }
    file_.print("ITEM: TIMESTEP\n{:d}\n", frame.step());
    file_.print("ITEM: NUMBER OF ATOMS\n{:d}\n", frame.size());

    const auto& cell = frame.cell();
    if (cell.shape() == UnitCell::ORTHORHOMBIC || cell.shape() == UnitCell::INFINITE) {
        file_.print("ITEM: BOX BOUNDS pp pp pp\n");
        auto lengths = cell.lengths();
        // print zeros if the cell is infinite, this line is still required
        file_.print("{:-1.12e} {:-1.12e}\n", 0.0, lengths[0]);
        file_.print("{:-1.12e} {:-1.12e}\n", 0.0, lengths[1]);
        file_.print("{:-1.12e} {:-1.12e}\n", 0.0, lengths[2]);
    } else { // Triclinic
        const auto& matrix = cell.matrix();
        if (!is_upper_triangular(matrix)) {
            throw format_error("unsupported triclinic but non upper-triangular cell "
                               "matrix in LAMMPS writer");
        }
        file_.print("ITEM: BOX BOUNDS xy xz yz pp pp pp\n");
        file_.print("{:-1.12e} {:-1.12e} {:-1.12e}\n", 0.0, matrix[0][0], matrix[0][1]);
        file_.print("{:-1.12e} {:-1.12e} {:-1.12e}\n", 0.0, matrix[1][1], matrix[0][2]);
        file_.print("{:-1.12e} {:-1.12e} {:-1.12e}\n", 0.0, matrix[2][2], matrix[1][2]);
    }

    bool has_names = false;
    for (size_t i = 0; i < frame.size(); ++i) {
        const auto& atom = frame[i];
        if (!atom.name().empty()) {
            has_names = true;
            break;
        }
    }

    auto positions = frame.positions();
    auto velocities = frame.velocities();
    file_.print("ITEM: ATOMS id xu yu zu type"); // write unwrapped positions
    if (has_names) {
        file_.print(" element");
    }
    file_.print(" mass q");
    if (velocities) {
        file_.print(" vx vy vz");
    }
    file_.print("\n");
    for (size_t i = 0; i < frame.size(); ++i) {
        const auto& atom = frame[i];
        // LAMMPS uses atom IDs that start with 1
        file_.print("{:d} {:g} {:g} {:g}", i + 1, positions[i][0], positions[i][1],
                    positions[i][2]);
        auto type = parse_lammps_type(atom.type());
        if (type && (min_numeric_type_ == 0 || *type <= min_numeric_type_)) {
            // a valid numeric type and no other invalid types encountered previously
            file_.print(" {:d}", *type);
            if (*type > max_numeric_type_) {
                max_numeric_type_ = *type;
            }
        } else {
            // use a generated numeric type because trajectory contains invalid types
            auto search = type_list_.find(atom.type());
            if (search != type_list_.end()) {
                // type has already a numeric type
                file_.print(" {:d}", search->second);
            } else {
                // a new type to generate a numeric type for
                min_numeric_type_ = max_numeric_type_;
                ++max_numeric_type_;
                type_list_.emplace(atom.type(), max_numeric_type_);
                file_.print(" {:d}", max_numeric_type_);
                warning("LAMMPS writer",
                        "trajectory with invalid types: generated type for '{}' is {}", atom.type(),
                        max_numeric_type_);
            }
        }
        if (has_names) {
            file_.print(" {:s}", atom.name());
        }
        file_.print(" {:g} {:g}", atom.mass(), atom.charge());
        if (velocities) {
            const auto& v = (*velocities)[i];
            file_.print(" {:g} {:g} {:g}", v[0], v[1], v[2]);
        }
        file_.print("\n");
    }
}

optional<uint64_t> LAMMPSTrajectoryFormat::forward() {
    auto position = file_.tellpos();
    size_t natoms = 0;

    auto line = file_.readline();
    if (trim(line).empty() || file_.eof()) {
        // no more data to read, so give up here
        return nullopt;
    }
    auto item = get_item(line);
    while (!file_.eof() && (!item || *item != "NUMBER OF ATOMS")) {
        line = file_.readline();
        item = get_item(line);
    }
    if (!item || *item != "NUMBER OF ATOMS") {
        throw format_error(
            "could not find 'ITEM: NUMBER OF ATOMS' in LAMMPS trajectory"
        );
    }

    try {
        line = file_.readline();
        natoms = parse<size_t>(trim(line));
    } catch (const Error&) {
        throw format_error(
            "could not parse the number of atoms in '{}' for LAMMPS trajectory",
            line
        );
    }

    // read the box
    for (size_t i = 0; i < 4; ++i) {
        line = file_.readline();
        if (file_.eof()) {
            throw format_error(
                "could not read box data in LAMMPS trajectory: not enough lines in the file"
            );
        }

        if (i == 0) {
            item = get_item(line);
            if (!item || item->substr(0, 10) != "BOX BOUNDS") {
                throw format_error(
                    "expected 'ITEM: BOX BOUNDS' after the number of atoms in "
                    "LAMMPS trajectory, got '{}'", line
                );
            }
        }
    }

    line = file_.readline();
    item = get_item(line);
    if (!item || item->substr(0, 5) != "ATOMS") {
        throw format_error(
            "could not read atom header for LAMMPS trajectory in this line: '{}'",
            line
        );
    }

    for (size_t i = 0; i < natoms; ++i) {
        file_.readline();
        if (file_.eof()) {
            throw format_error(
                "this file does not contain enough lines in ATOMS section for LAMMPS trajectory"
            );
        }
    }

    return position;
}
