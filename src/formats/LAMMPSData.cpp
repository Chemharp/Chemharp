// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cstddef>
#include <cstdint>
#include <cassert>
#include <cmath>

#include <tuple>
#include <array>
#include <string>
#include <vector>
#include <algorithm>
#include <string_view>
#include <unordered_map>
#include <unordered_set>

#include "chemfiles/types.hpp"
#include "chemfiles/utils.hpp"
#include "chemfiles/parse.hpp"
#include "chemfiles/warnings.hpp"
#include "chemfiles/error_fmt.hpp"
#include "chemfiles/sorted_set.hpp"
#include "chemfiles/unreachable.hpp"
#include "chemfiles/external/span.hpp"
#include "chemfiles/external/optional.hpp"

#include "chemfiles/Atom.hpp"
#include "chemfiles/File.hpp"
#include "chemfiles/Frame.hpp"
#include "chemfiles/Residue.hpp"
#include "chemfiles/Topology.hpp"
#include "chemfiles/UnitCell.hpp"
#include "chemfiles/Connectivity.hpp"
#include "chemfiles/FormatMetadata.hpp"

#include "chemfiles/formats/LAMMPSData.hpp"


using namespace chemfiles;

/// Guess molecules id from the bonds. This function return a vector containing
/// the molecule id for each atom in the frame
static std::vector<size_t> guess_molecules(const Frame& frame);

/// Make sure the tilt factor matrix[i][j] is contained between -matrix[i][i] / 2
/// and matrix[i][i] / 2.
static double tilt_factor(const Matrix3D& matrix, size_t i, size_t j);

template<> const FormatMetadata& chemfiles::format_metadata<LAMMPSDataFormat>() {
    static FormatMetadata metadata;
    metadata.name = "LAMMPS Data";
    metadata.extension = nullopt;
    metadata.description = "LAMMPS input data file";
    metadata.reference = "https://lammps.sandia.gov/doc/read_data.html";

    metadata.read = true;
    metadata.write = true;
    metadata.memory = true;

    metadata.positions = true;
    metadata.velocities = true;
    metadata.unit_cell = true;
    metadata.atoms = true;
    metadata.bonds = true;
    metadata.residues = true;
    return metadata;
}

atom_style::atom_style(std::string name): name_(std::move(name)) {
    if (name_ == "angle") {
        style_ = ANGLE;
    } else if (name_ == "atomic") {
        style_ = ATOMIC;
    } else if (name_ == "body") {
        style_ = BODY;
    } else if (name_ == "bond") {
        style_ = BOND;
    } else if (name_ == "charge") {
        style_ = CHARGE;
    } else if (name_ == "dipole") {
        style_ = DIPOLE;
    } else if (name_ == "dpd") {
        style_ = DPD;
    } else if (name_ == "electron") {
        style_ = ELECTRON;
    } else if (name_ == "ellipsoid") {
        style_ = ELLIPSOID;
    } else if (name_ == "full") {
        style_ = FULL;
    } else if (name_ == "line") {
        style_ = LINE;
    } else if (name_ == "meso") {
        style_ = MESO;
    } else if (name_ == "molecular") {
        style_ = MOLECULAR;
    } else if (name_ == "peri") {
        style_ = PERI;
    } else if (name_ == "smd") {
        style_ = SMD;
    } else if (name_ == "sphere") {
        style_ = SPHERE;
    } else if (name_ == "template") {
        style_ = TEMPLATE;
    } else if (name_ == "tri") {
        style_ = TRI;
    } else if (name_ == "wavepacket") {
        style_ = WAVEPACKET;
    } else if (name_ == "hybrid") {
        style_ = HYBRID;
    } else {
        throw format_error("LAMMPS Data: unknown atom style '{}'", name_);
    }
}

atom_data atom_style::read_line(std::string_view line, size_t index) const {
    atom_data d;

    // dummy variable to ignore a value
    double SKIP;
    switch (style_) {
    case ANGLE:
    case BOND:
    case MOLECULAR:
        // atom-ID molecule-ID atom-type x y z
        scan(line, d.index, d.molid, d.type, d.x, d.y, d.z);
        break;
    case ATOMIC:
        // atom-ID atom-type x y z
        scan(line, d.index, d.type, d.x, d.y, d.z);
        break;
    case BODY:
        // atom-ID atom-type bodyflag mass x y z
        scan(line, d.index, d.type, SKIP, d.mass, d.x, d.y, d.z);
        break;
    case CHARGE:
    case DIPOLE:
        // atom-ID atom-type q x y z
        // atom-ID atom-type q x y z mux muy muz
        scan(line, d.index, d.type, d.charge, d.x, d.y, d.z);
        break;
    case DPD:
        // atom-ID atom-type theta x y z
        scan(line, d.index, d.type, SKIP, d.x, d.y, d.z);
        break;
    case ELECTRON:
        // atom-ID atom-type q spin eradius x y z
        scan(line, d.index, d.type, SKIP, SKIP, d.x, d.y, d.z);
        break;
    case ELLIPSOID:
        // atom-ID atom-type ellipsoidflag density x y z
        scan(line, d.index, d.type, SKIP, SKIP, d.x, d.y, d.z);
        break;
    case FULL:
        // atom-ID molecule-ID atom-type q x y z
        scan(line, d.index, d.molid, d.type, d.charge, d.x, d.y, d.z);
        break;
    case LINE:
        // atom-ID molecule-ID atom-type lineflag density x y z
        scan(line, d.index, d.molid, d.type, SKIP, SKIP, d.x, d.y, d.z);
        break;
    case MESO:
        // atom-ID atom-type rho e cv x y z
        scan(line, d.index, d.type, SKIP, SKIP, SKIP, d.x, d.y, d.z);
        break;
    case PERI:
        // atom-ID atom-type volume density x y z
        scan(line, d.index, d.type, SKIP, SKIP, d.x, d.y, d.z);
        break;
    case SMD:
        // atom-ID atom-type molecule volume mass kernel-radius contact-radius x y z
        scan(line, d.index, d.type, d.molid, SKIP, d.mass, SKIP, SKIP, d.x, d.y, d.z);
        break;
    case SPHERE:
        // atom-ID atom-type diameter density x y z
        scan(line, d.index, d.type, SKIP, SKIP, d.x, d.y, d.z);
        break;
    case TEMPLATE:
        // atom-ID molecule-ID template-index template-atom atom-type x y z
        scan(line, d.index, d.molid, SKIP, SKIP, d.type, d.x, d.y, d.z);
        break;
    case TRI:
        // atom-ID molecule-ID atom-type triangleflag density x y z
        scan(line, d.index, d.molid, d.type, SKIP, SKIP, d.x, d.y, d.z);
        break;
    case WAVEPACKET:
        // atom-ID atom-type charge spin eradius etag cs_re cs_im x y z
        scan(line, d.index, d.type, d.charge, SKIP, SKIP, SKIP, SKIP, SKIP, d.x, d.y, d.z);
        break;
    case HYBRID:
        if (!warned_) {
            warning("LAMMPS Data reader", "only reading the first style for atom_style hybrid");
            warned_ = true;
        }
        // atom-ID atom-type x y z sub-style1 sub-style2 ...
        scan(line, d.index, d.type, d.x, d.y, d.z);
        break;
    }

    if (d.index == 0) {
        // 0 means the user do not care about indexes, but we still need one.
        // So we use the index provided by the caller of this function
        d.index = index;
    } else {
        // LAMMPS uses 1 based indexes, convert it to 0-based
        d.index -= 1;
    }

    return d;
}

/// Remove the comment from `line` and return it.
static std::string_view split_comment(std::string_view& line);
/// Check if the line is an unused header value
static bool is_unused_header(std::string_view line);

void LAMMPSDataFormat::read_next(Frame& frame) {
    if (file_.tellpos() != 0) {
        throw format_error("LAMMPS Data format only supports reading one frame");
    }

    auto comment = file_.readline();
    // VMD topotools writes the atom style in the comment header
    auto it = comment.find("atom_style");
    if (it != std::string::npos) {
        auto style = comment.substr(it + 10);
        atom_style_name_ = std::string(trim(split(style, ' ')[0]));
    }

    while(!file_.eof()) {
        switch (current_section_) {
        case HEADER:
            read_header(frame);
            break;
        case ATOMS:
            read_atoms(frame);
            break;
        case MASSES:
            read_masses();
            break;
        case BONDS:
            read_bonds(frame);
            break;
        case ANGLES:
            read_angles(frame);
            break;
        case DIHEDRALS:
            read_dihedrals(frame);
            break;
        case IMPROPERS:
            read_impropers(frame);
            break;
        case VELOCITIES:
            read_velocities(frame);
            break;
        case IGNORED:
            skip_to_next_section();
            break;
        case NOT_A_SECTION:
            unreachable();
        }
    }

    setup_masses(frame);
    setup_names(frame);
}

void LAMMPSDataFormat::read_header(Frame& frame) {
    assert(current_section_ == HEADER);
    auto matrix = Matrix3D::unit();
    auto shape = UnitCell::ORTHORHOMBIC;

    while (!file_.eof()) {
        auto line = file_.readline();
        auto content = line;
        split_comment(content);

        if (content.empty() || is_unused_header(content)) {
            // Nothing to do
        } else if (content.find("atoms") != std::string::npos) {
            natoms_ = read_header_integer(content, "atoms");
        } else if (content.find("bonds") != std::string::npos) {
            nbonds_ = read_header_integer(content, "bonds");
        } else if (content.find("angles") != std::string::npos) {
            nangles_ = read_header_integer(content, "angles");
        } else if (content.find("dihedrals") != std::string::npos) {
            ndihedrals_ = read_header_integer(content, "dihedrals");
        } else if (content.find("impropers") != std::string::npos) {
            nimpropers_ = read_header_integer(content, "impropers");
        } else if (content.find("atom types") != std::string::npos) {
            natom_types_ = read_header_integer(content, "atom types");
        } else if (content.find("xlo") != std::string::npos) {
            matrix[0][0] = read_header_box_bounds(content, "xlo", "xhi");
        } else if (content.find("ylo") != std::string::npos) {
            matrix[1][1] = read_header_box_bounds(content, "ylo", "yhi");
        } else if (content.find("zlo") != std::string::npos) {
            matrix[2][2] = read_header_box_bounds(content, "zlo", "zhi");
        } else if (content.find("xy") != std::string::npos) {
            auto splitted = split(content, ' ');
            if (splitted.size() != 6 || splitted[3] != "xy" || splitted[4] != "xz" || splitted[5] != "yz") {
                throw format_error(
                    "invalid header value: expected '<xy> <xz> <yz> xy xz yz', got '{}'", content
                );
            }
            matrix[0][1] = parse<double>(splitted[0]);
            matrix[0][2] = parse<double>(splitted[1]);
            matrix[1][2] = parse<double>(splitted[2]);
            // Even if all parameters are 0, set shape to TRICLINIC
            shape = UnitCell::TRICLINIC;
        } else {
            // End of the header, get the section and break
            current_section_ = get_section(line);
            assert(current_section_ != NOT_A_SECTION);
            break;
        }
    }
    auto cell = UnitCell(matrix);
    cell.set_shape(shape);
    frame.set_cell(cell);
}

size_t LAMMPSDataFormat::read_header_integer(std::string_view line, const std::string& context) {
    auto splitted = split(line, ' ');
    if (splitted.size() < 2) {
        throw format_error(
            "invalid header value: expected '<n> {}', got '{}'", context, line
        );
    }
    return parse<size_t>(splitted[0]);
}

double LAMMPSDataFormat::read_header_box_bounds(std::string_view line, const std::string& lo, const std::string& hi) {
    auto splitted = split(line, ' ');
    if (splitted.size() < 4 || splitted[2] != lo || splitted[3] != hi) {
        throw format_error(
            "invalid header value: expected '<lo> <hi> {} {}', got '{}'", lo, hi, line
        );
    }
    auto low = parse<double>(splitted[0]);
    auto high = parse<double>(splitted[1]);
    return high - low;
}

void LAMMPSDataFormat::get_next_section() {
    while (!file_.eof()) {
        auto line = file_.readline();
        if (!line.empty()) {
            auto section = get_section(line);
            if (section == NOT_A_SECTION) {
                throw format_error("expected section name, got '{}'", line);
            } else {
                current_section_ = section;
                break;
            }
        }
    }
}

void LAMMPSDataFormat::skip_to_next_section() {
    while (!file_.eof()) {
        auto line = file_.readline();
        if (!line.empty()) {
            auto section = get_section(line);
            if (section == NOT_A_SECTION) {
                continue;
            } else {
                current_section_ = section;
                break;
            }
        }
    }
}

void LAMMPSDataFormat::read_atoms(Frame& frame) {
    assert(current_section_ == ATOMS);
    if (natoms_ == 0) {
        throw format_error("missing atoms count in header");
    }

    if (atom_style_name_.empty()) {
        warning("LAMMPS Data reader", "unknown atom style, defaulting to 'full'");
        atom_style_name_ = "full";
    }
    style_ = atom_style(atom_style_name_);

    frame.resize(natoms_);
    auto positions = frame.positions();
    auto residues = std::unordered_map<size_t, Residue>();

    size_t n = 0;
    while (n < natoms_ && !file_.eof()) {
        auto line = file_.readline();
        auto comment = split_comment(line);
        if (line.empty()) {
            continue;
        }

        auto data = style_.read_line(line, n);
        if (data.index >= natoms_) {
            throw format_error(
                "too many atoms in [Atoms] section: expected {} atoms, got atom with index {}",
                natoms_, data.index
            );
        }

        if (!comment.empty()) {
            // Read the first string after the comment, and use it as atom name
            auto name = split(comment, ' ')[0];
            if (names_.empty()) {
                names_.resize(natoms_);
            }
            names_[data.index] = std::string(name);
        }

        auto atom = Atom(std::to_string(data.type));
        if (!std::isnan(data.charge)) {
            atom.set_charge(data.charge);
        }
        if (!std::isnan(data.mass)) {
            atom.set_mass(data.mass);
        }

        if (data.molid != 0) {
            auto residue_iter = residues.find(data.molid);
            if (residue_iter != residues.end()) {
                residue_iter->second.add_atom(data.index);
            } else {
                auto residue = Residue("", static_cast<int64_t>(data.molid));
                residue.add_atom(data.index);
                residues.emplace(data.molid, std::move(residue));
            }
        }

        frame[data.index] = atom;
        positions[data.index][0] = data.x;
        positions[data.index][1] = data.y;
        positions[data.index][2] = data.z;
        n++;
    }

    for (auto it: std::move(residues)) {
        frame.add_residue(std::move(it.second));
    }

    get_next_section();
}

void LAMMPSDataFormat::read_masses() {
    assert(current_section_ == MASSES);
    if (natom_types_ == 0) {
        throw format_error("missing atom types count in header");
    }
    size_t n = 0;
    while (n < natom_types_ && !file_.eof()) {
        auto line = file_.readline();
        split_comment(line);
        if (line.empty()) {
            continue;
        }

        auto splitted = split(line, ' ');
        if (splitted.size() != 2) {
            throw format_error("bad mass specification '{}'", line);
        }

        auto type = splitted[0];
        auto mass = parse<double>(splitted[1]);
        masses_.emplace(std::string(type), mass);
        n++;
    }

    get_next_section();
}

void LAMMPSDataFormat::read_bonds(Frame& frame) {
    assert(current_section_ == BONDS);
    if (nbonds_ == 0) {
        throw format_error("missing bonds count in header");
    }
    size_t n = 0;
    while (n < nbonds_ && !file_.eof()) {
        auto line = file_.readline();
        split_comment(line);
        if (line.empty()) {
            continue;
        }

        auto splitted = split(line, ' ');
        if (splitted.size() != 4) {
            throw format_error("bad bond specification '{}'", line);
        }
        // LAMMPS use 1-based indexing
        auto bond_type = parse<std::string>(splitted[1]);
        auto i = parse<size_t>(splitted[2]) - 1;
        auto j = parse<size_t>(splitted[3]) - 1;
        frame.add_bond(i, j, Bond::UNKNOWN, bond_type);
        n++;
    }

    if (file_.eof() && n < nbonds_) {
        throw format_error("end of file found before getting all bonds");
    }

    get_next_section();
}

void LAMMPSDataFormat::read_angles(Frame& frame) {
    assert(current_section_ == ANGLES);
    if (nangles_ == 0) {
        throw format_error("missing angles count in header");
    }
    size_t n = 0;
    while (n < nangles_ && !file_.eof()) {
        auto line = file_.readline();
        split_comment(line);
        if (line.empty()) {continue;}

        auto splitted = split(line, ' ');
        if (splitted.size() != 5) {
            throw format_error("bad angle specification '{}'", line);
        }
        // LAMMPS use 1-based indexing
        auto angle_type = parse<std::string>(splitted[1]);
        auto i = parse<size_t>(splitted[2]) - 1;
        auto j = parse<size_t>(splitted[3]) - 1;
        auto k = parse<size_t>(splitted[4]) - 1;
        frame.add_angle(i, j, k, angle_type);
        n++;
    }

    if (file_.eof() && n < nangles_) {
        throw format_error("end of file found before getting all angles");
    }

    get_next_section();
}

void LAMMPSDataFormat::read_dihedrals(Frame& frame) {
    assert(current_section_ == DIHEDRALS);
    if (ndihedrals_ == 0) {
        throw format_error("missing dihedrals count in header");
    }
    size_t n = 0;
    while (n < ndihedrals_ && !file_.eof()) {
        auto line = file_.readline();
        split_comment(line);
        if (line.empty()) {continue;}

        auto splitted = split(line, ' ');
        if (splitted.size() != 6) {
            throw format_error("bad dihedral specification '{}'", line);
        }
        // LAMMPS use 1-based indexing
        auto dihedral_type = parse<std::string>(splitted[1]);
        auto i = parse<size_t>(splitted[2]) - 1;
        auto j = parse<size_t>(splitted[3]) - 1;
        auto k = parse<size_t>(splitted[4]) - 1;
        auto m = parse<size_t>(splitted[5]) - 1;
        frame.add_dihedral(i, j, k, m, dihedral_type);
        n++;
    }

    if (file_.eof() && n < ndihedrals_) {
        throw format_error("end of file found before getting all dihedrals");
    }

    get_next_section();
}

void LAMMPSDataFormat::read_impropers(Frame& frame) {
    assert(current_section_ == IMPROPERS);
    if (nimpropers_ == 0) {
        throw format_error("missing impropers count in header");
    }
    size_t n = 0;
    while (n < nimpropers_ && !file_.eof()) {
        auto line = file_.readline();
        split_comment(line);
        if (line.empty()) {continue;}

        auto splitted = split(line, ' ');
        if (splitted.size() != 6) {
            throw format_error("bad improper specification '{}'", line);
        }
        // LAMMPS use 1-based indexing
        auto improper_type = parse<std::string>(splitted[1]);
        auto i = parse<size_t>(splitted[2]) - 1;
        auto j = parse<size_t>(splitted[3]) - 1;
        auto k = parse<size_t>(splitted[4]) - 1;
        auto m = parse<size_t>(splitted[5]) - 1;
        frame.add_improper(i, j, k, m, improper_type);
        n++;
    }

    if (file_.eof() && n < nimpropers_) {
        throw format_error("end of file found before getting all impropers");
    }

    get_next_section();
}

void LAMMPSDataFormat::read_velocities(Frame& frame) {
    assert(current_section_ == VELOCITIES);
    if (natoms_ == 0) {
        throw format_error("missing atoms count in header");
    }
    size_t n = 0;
    frame.add_velocities();
    auto velocities = *frame.velocities();
    while (n < natoms_ && !file_.eof()) {
        auto line = file_.readline();
        split_comment(line);
        if (line.empty()) {
            continue;
        }

        auto splitted = split(line, ' ');
        if (splitted.size() < 4) {
            throw format_error("bad velocity specification '{}'", line);
        }
        // LAMMPS use 1-based indexing
        auto id = parse<size_t>(splitted[0]) - 1;
        auto vx = parse<double>(splitted[1]);
        auto vy = parse<double>(splitted[2]);
        auto vz = parse<double>(splitted[3]);
        velocities[id] = Vector3D(vx, vy, vz);
        n++;
    }

    if (file_.eof() && n < nbonds_) {
        throw format_error("end of file found before getting all velocities");
    }

    get_next_section();
}


void LAMMPSDataFormat::setup_masses(Frame& frame) const {
    if (masses_.empty()) {return;}
    for (auto& atom: frame) {
        auto it = masses_.find(atom.type());
        if (it != masses_.end()) {
            atom.set_mass(it->second);
        }
    }
}

void LAMMPSDataFormat::setup_names(Frame& frame) const {
    if (names_.empty()) {return;}
    assert(names_.size() == frame.size());

    for (size_t i=0; i<frame.size(); i++) {
        if (!names_[i].empty()) {
            frame[i].set_name(names_[i]);
            frame[i].set_type(names_[i]);
        }
    }
}


static std::unordered_set<std::string_view> IGNORED_SECTIONS = {
    "Ellipsoids", "Lines", "Triangles", "Bodies", "Pair Coeffs",
    "PairIJ Coeffs", "Bond Coeffs", "Angle Coeffs", "Dihedral Coeffs",
    "Improper Coeffs", "BondBond Coeffs", "BondAngle Coeffs",
    "MiddleBondTorsion Coeffs", "EndBondTorsion Coeffs", "AngleTorsion Coeffs",
    "AngleAngleTorsion Coeffs", "BondBond13 Coeffs", "AngleAngle Coeffs"
};

LAMMPSDataFormat::section_t LAMMPSDataFormat::get_section(std::string_view section) {
    auto comment = split_comment(section);
    if (section == "Atoms") {
        if (!comment.empty()) {
            atom_style_name_ = std::string(trim(comment));
        }
        return ATOMS;
    } else if (section == "Bonds") {
        return BONDS;
    } else if (section == "Velocities") {
        return VELOCITIES;
    } else if (section == "Masses") {
        return MASSES;
    } else if (section == "Angles" || section == "Dihedrals" || section == "Impropers") {
        // We don't use the angles, dihedral and improper sections, but we don't
        // send a warning for them, as they are guessed from the bonds

        // TODO: check that the angles and dihedral list matches
        return IGNORED;
    } else if (IGNORED_SECTIONS.find(section) != IGNORED_SECTIONS.end()) {
        if (section.find("Coeffs") == std::string::npos) {
            // Don't send a warning for force field parameters
            warning("LAMMPS Data reader", "ignoring section '{}'", section);
        }
        return IGNORED;
    } else {
        return NOT_A_SECTION;
    }
}

static bond_type normalize_bond_type(size_t i, size_t j) {
    if (i < j) {
        return std::make_tuple(i, j);
    } else {
        return std::make_tuple(j, i);
    }
}

static angle_type normalize_angle_type(size_t i, size_t j, size_t k) {
    if (i < k) {
        return std::make_tuple(i, j, k);
    } else {
        return std::make_tuple(k, j, i);
    }
}

static dihedral_type normalize_dihedral_type(size_t i, size_t j, size_t k, size_t m) {
    auto max_ij = std::max(i, j);
    auto max_km = std::max(k, m);
    if (max_ij == max_km) {
        if (std::min(i, j) < std::min(k, m)) {
            return std::make_tuple(i, j, k, m);
        } else {
            return std::make_tuple(m, k, j, i);
        }
    } else if (max_ij < max_km) {
        return std::make_tuple(i, j, k, m);
    } else {
        return std::make_tuple(m, k, j, i);
    }
}

static improper_type normalize_improper_type(size_t i, size_t j, size_t k, size_t m) {
    std::array<size_t, 3> others = {{i, k, m}};
    std::sort(others.begin(), others.end());
    i = others[0];
    k = others[1];
    m = others[2];
    return std::make_tuple(i, j, k, m);
}

static atom_type make_atom_type(const Atom& atom) {
    return {atom.type(), atom.mass()};
}

DataTypes::DataTypes(const Topology& topology) {
    for (auto& atom: topology) {
        atoms_.insert(make_atom_type(atom));
    }

    for (auto& bond: topology.bonds()) {
        auto i = atom_type_id(topology[bond[0]]);
        auto j = atom_type_id(topology[bond[1]]);
        bonds_.insert({normalize_bond_type(i, j), topology.bond_type(bond[0], bond[1])});
    }

    for (auto& angle: topology.angles()) {
        auto i = atom_type_id(topology[angle[0]]);
        auto j = atom_type_id(topology[angle[1]]);
        auto k = atom_type_id(topology[angle[2]]);
        angles_.insert({normalize_angle_type(i, j, k), topology.angle_type(angle[0], angle[1], angle[2])});
    }

    for (auto& dihedral: topology.dihedrals()) {
        auto i = atom_type_id(topology[dihedral[0]]);
        auto j = atom_type_id(topology[dihedral[1]]);
        auto k = atom_type_id(topology[dihedral[2]]);
        auto m = atom_type_id(topology[dihedral[3]]);
        dihedrals_.insert({normalize_dihedral_type(i, j, k, m), topology.dihedral_type(dihedral[0], dihedral[1], dihedral[2], dihedral[3])});
    }

    for (auto& improper: topology.impropers()) {
        auto i = atom_type_id(topology[improper[0]]);
        auto j = atom_type_id(topology[improper[1]]);
        auto k = atom_type_id(topology[improper[2]]);
        auto m = atom_type_id(topology[improper[3]]);
        impropers_.insert({normalize_improper_type(i, j, k, m), topology.improper_type(improper[0], improper[1], improper[2], improper[3])});
    }
}

size_t DataTypes::atom_type_id(const Atom& atom) const {
    auto it = atoms_.find(make_atom_type(atom));
    if (it != atoms_.end()) {
        return static_cast<size_t>(it - atoms_.begin());
    } else {
        throw error("invalid atom type passed to atom_type_id, this is a bug");
    }
}

std::string DataTypes::bond_type_id(size_t type_i, size_t type_j) const {
    auto it = std::find_if(bonds_.begin(), bonds_.end(), [&type_i, &type_j](const std::tuple<bond_type, std::string>& pair) {
        return std::get<0>(pair) == normalize_bond_type(type_i, type_j);
    });
    if (it != bonds_.end()) {
        if (std::get<1>(*it).empty()) {
            return std::to_string(it - bonds_.begin() + 1);
        }
        return std::get<1>(*it);
    } else {
        throw error("invalid bond type passed to bond_type_id, this is a bug");
    }
}

std::string DataTypes::angle_type_id(size_t type_i, size_t type_j, size_t type_k) const {
    auto it = std::find_if(angles_.begin(), angles_.end(), [&type_i, &type_j, &type_k](const std::tuple<angle_type, std::string>& pair) {
        return std::get<0>(pair) == normalize_angle_type(type_i, type_j, type_k);
    });
    if (it != angles_.end()) {
        if (std::get<1>(*it).empty()) {
            return std::to_string(it - angles_.begin() + 1);
        }
        return std::get<1>(*it);
    } else {
        throw error("invalid angle type passed to angle_type_id, this is a bug");
    }
}

std::string DataTypes::dihedral_type_id(size_t type_i, size_t type_j, size_t type_k, size_t type_m) const {
    auto it = std::find_if(dihedrals_.begin(), dihedrals_.end(), [&type_i, &type_j, &type_k, &type_m](const std::tuple<dihedral_type, std::string>& pair) {
        return std::get<0>(pair) == normalize_dihedral_type(type_i, type_j, type_k, type_m);
    });
    if (it != dihedrals_.end()) {
        if (std::get<1>(*it).empty()) {
            return std::to_string(it - dihedrals_.begin() + 1);
        }
        return std::get<1>(*it);
    } else {
        throw error("invalid dihedral type passed to dihedral_type_id, this is a bug");
    }
}

std::string DataTypes::improper_type_id(size_t type_i, size_t type_j, size_t type_k, size_t type_m) const {
    auto it = std::find_if(impropers_.begin(), impropers_.end(), [&type_i, &type_j, &type_k, &type_m](const std::tuple<improper_type, std::string>& pair) {
        return std::get<0>(pair) == normalize_improper_type(type_i, type_j, type_k, type_m);
    });
    if (it != impropers_.end()) {
        if (std::get<1>(*it).empty()) {
            return std::to_string(it - impropers_.begin() + 1);
        }
        return std::get<1>(*it);
    } else {
        throw error("invalid improper type passed to improper_type_id, this is a bug");
    }
}

void LAMMPSDataFormat::write_next(const Frame& frame) {
    if (file_.tellpos() != 0) {
        throw format_error("LAMMPS Data format only supports writing one frame");
    }

    auto types = DataTypes(frame.topology());

    write_header(types, frame);
    write_types(types);
    write_masses(types);
    write_atoms(types, frame);
    write_velocities(frame);
    auto& topology = frame.topology();
    write_bonds(types, topology);
    write_angles(types, topology);
    write_dihedrals(types, topology);
    write_impropers(types, topology);
}

void LAMMPSDataFormat::write_header(const DataTypes& types, const Frame& frame) {
    file_.print("LAMMPS data file -- atom_style full -- generated by chemfiles\n\n");
    file_.print("{} atoms\n", frame.size());
    file_.print("{} bonds\n", frame.topology().bonds().size());
    file_.print("{} angles\n", frame.topology().angles().size());
    file_.print("{} dihedrals\n", frame.topology().dihedrals().size());
    file_.print("{} impropers\n", frame.topology().impropers().size());

    file_.print("{} atom types\n", types.atoms().size());
    file_.print("{} bond types\n", types.bonds().size());
    file_.print("{} angle types\n", types.angles().size());
    file_.print("{} dihedral types\n", types.dihedrals().size());
    file_.print("{} improper types\n", types.impropers().size());

    auto matrix = frame.cell().matrix();
    file_.print("{:#.9} {:#.9} xlo xhi\n", 0.0, matrix[0][0]);
    file_.print("{:#.9} {:#.9} ylo yhi\n", 0.0, matrix[1][1]);
    file_.print("{:#.9} {:#.9} zlo zhi\n", 0.0, matrix[2][2]);
    if (frame.cell().shape() == UnitCell::TRICLINIC) {
        assert(tilt_factor(matrix, 1, 0) == 0);
        assert(tilt_factor(matrix, 2, 0) == 0);
        assert(tilt_factor(matrix, 2, 1) == 0);
        file_.print("{:#.9} {:#.9} {:#.9} xy xz yz\n",
            tilt_factor(matrix, 0, 1),
            tilt_factor(matrix, 0, 2),
            tilt_factor(matrix, 1, 2)
        );
    }

    file_.print("\n");
}

void LAMMPSDataFormat::write_types(const DataTypes& types) {
    auto& atoms = types.atoms().as_vec();
    if (!atoms.empty()) {
        file_.print("# Pair Coeffs\n");
        for (size_t i=0; i<atoms.size(); i++) {
            file_.print("# {} {}\n", i + 1, atoms[i].first);
        }
    }

    auto& bonds = types.bonds().as_vec();
    if (!bonds.empty()) {
        file_.print("\n# Bond Coeffs\n");
        for (size_t i=0; i<bonds.size(); i++) {
            auto& bond_index = std::get<0>(bonds[i]);
            file_.print("# {} {}-{}\n", i + 1,
                atoms[std::get<0>(bond_index)].first,
                atoms[std::get<1>(bond_index)].first
            );
        }
    }

    auto& angles = types.angles().as_vec();
    if (!angles.empty()) {
        file_.print("\n# Angle Coeffs\n");
        for (size_t i=0; i<angles.size(); i++) {
            auto& angle_index = std::get<0>(angles[i]);
            file_.print("# {} {}-{}-{}\n", i + 1,
                atoms[std::get<0>(angle_index)].first,
                atoms[std::get<1>(angle_index)].first,
                atoms[std::get<2>(angle_index)].first
            );
        }
    }

    auto& dihedrals = types.dihedrals().as_vec();
    if (!dihedrals.empty()) {
        file_.print("\n# Dihedrals Coeffs\n");
        for (size_t i=0; i<dihedrals.size(); i++) {
            auto& dihedral_index = std::get<0>(dihedrals[i]);
            file_.print("# {} {}-{}-{}-{}\n", i + 1,
                atoms[std::get<0>(dihedral_index)].first,
                atoms[std::get<1>(dihedral_index)].first,
                atoms[std::get<2>(dihedral_index)].first,
                atoms[std::get<3>(dihedral_index)].first
            );
        }
    }

    auto& impropers = types.impropers().as_vec();
    if (!impropers.empty()) {
        file_.print("\n# Impropers Coeffs\n");
        for (size_t i=0; i<impropers.size(); i++) {
            auto& improper_index = std::get<0>(impropers[i]);
            file_.print("# {} {}-{}-{}-{}\n", i + 1,
                atoms[std::get<0>(improper_index)].first,
                atoms[std::get<1>(improper_index)].first,
                atoms[std::get<2>(improper_index)].first,
                atoms[std::get<3>(improper_index)].first
            );
        }
    }
}

void LAMMPSDataFormat::write_masses(const DataTypes& types) {
    file_.print("\nMasses\n\n");
    auto& atoms = types.atoms().as_vec();
    for (size_t i=0; i<atoms.size(); i++) {
        file_.print("{} {:#g} # {}\n", i + 1, atoms[i].second, atoms[i].first);
    }
}

void LAMMPSDataFormat::write_atoms(const DataTypes& types, const Frame& frame) {
    file_.print("\nAtoms # full\n\n");
    auto positions = frame.positions();
    auto molids = guess_molecules(frame);
    for (size_t i=0; i<frame.size(); i++) {
        auto& atom = frame.topology()[i];
        auto molid = molids[i];
        file_.print("{} {} {} {:#g} {:#g} {:#g} {:#g} # {}\n",
            i + 1, molid + 1, types.atom_type_id(atom) + 1, atom.charge(),
            positions[i][0], positions[i][1], positions[i][2],
            atom.type()
        );
    }
}

void LAMMPSDataFormat::write_velocities(const Frame& frame) {
    if (!frame.velocities()) { return; }

    file_.print("\nVelocities\n\n");
    auto velocities = *frame.velocities();
    for (size_t i=0; i<frame.size(); i++) {
        file_.print("{} {} {} {}\n",
            i + 1, velocities[i][0], velocities[i][1], velocities[i][2]
        );
    }
}

void LAMMPSDataFormat::write_bonds(const DataTypes& types, const Topology& topology) {
    if (topology.bonds().empty()) { return; }

    file_.print("\nBonds\n\n");
    size_t bond_id = 1;
    for (auto bond: topology.bonds()) {
        auto type_i = types.atom_type_id(topology[bond[0]]);
        auto type_j = types.atom_type_id(topology[bond[1]]);
        auto bond_type_str = types.bond_type_id(type_i, type_j);
        size_t bond_type_id = parse<size_t>(bond_type_str);
        // if (bond_type_str.empty()) {
        //     bond_type_id += 1;
        // }
        file_.print("{} {} {} {}\n",
            bond_id, bond_type_id, bond[0] + 1, bond[1] + 1
        );
        bond_id++;
    }
}

void LAMMPSDataFormat::write_angles(const DataTypes& types, const Topology& topology) {
    if (topology.angles().empty()) { return; }

    file_.print("\nAngles\n\n");
    size_t angle_id = 1;
    for (auto angle: topology.angles()) {
        auto type_i = types.atom_type_id(topology[angle[0]]);
        auto type_j = types.atom_type_id(topology[angle[1]]);
        auto type_k = types.atom_type_id(topology[angle[2]]);
        auto angle_type_str = types.angle_type_id(type_i, type_j, type_k);
        size_t angle_type_id = parse<size_t>(angle_type_str);
        // if (angle_type_str.empty()) {
        //     angle_type_id += 1;
        // }
        file_.print("{} {} {} {} {}\n",
            angle_id, angle_type_id,
            angle[0] + 1, angle[1] + 1, angle[2] + 1
        );
        angle_id++;
    }
}

void LAMMPSDataFormat::write_dihedrals(const DataTypes& types, const Topology& topology) {
    if (topology.dihedrals().empty()) { return; }

    file_.print("\nDihedrals\n\n");
    size_t dihedral_id = 1;
    for (auto dihedral: topology.dihedrals()) {
        auto type_i = types.atom_type_id(topology[dihedral[0]]);
        auto type_j = types.atom_type_id(topology[dihedral[1]]);
        auto type_k = types.atom_type_id(topology[dihedral[2]]);
        auto type_m = types.atom_type_id(topology[dihedral[3]]);
        auto dihedral_type_str = types.dihedral_type_id(type_i, type_j, type_k, type_m);
        size_t dihedral_type_id = parse<size_t>(dihedral_type_str);
        // if (dihedral_type_str.empty()) {
        //     dihedral_type_id += 1;
        // }
        file_.print("{} {} {} {} {} {}\n",
            dihedral_id, dihedral_type_id,
            dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1
        );
        dihedral_id++;
    }
}

void LAMMPSDataFormat::write_impropers(const DataTypes& types, const Topology& topology) {
    if (topology.impropers().empty()) { return; }

    file_.print("\nImpropers\n\n");
    size_t improper_id = 1;
    for (auto improper: topology.impropers()) {
        auto type_i = types.atom_type_id(topology[improper[0]]);
        auto type_j = types.atom_type_id(topology[improper[1]]);
        auto type_k = types.atom_type_id(topology[improper[2]]);
        auto type_m = types.atom_type_id(topology[improper[3]]);
        auto improper_type_str = types.improper_type_id(type_i, type_j, type_k, type_m);
        size_t improper_type_id = parse<size_t>(improper_type_str);
        // if (improper_type_str.empty()) {
        //     improper_type_id += 1;
        // }
        file_.print("{} {} {} {} {} {}\n",
            improper_id, improper_type_id,
            improper[0] + 1, improper[1] + 1, improper[2] + 1, improper[3] + 1
        );
        improper_id++;
    }
}

std::string_view split_comment(std::string_view& line) {
    line = trim(line);
    auto position = line.find('#');
    if (position != std::string::npos) {
        auto comment = line.substr(position + 1);
        line.remove_suffix(line.size() - position);
        line = trim(line);
        return comment;
    } else {
        return "";
    }
}

bool is_unused_header(std::string_view line) {
    return (line.find("angles") != std::string::npos) ||
           (line.find("dihedrals") != std::string::npos) ||
           (line.find("impropers") != std::string::npos) ||
           (line.find("bond types") != std::string::npos) ||
           (line.find("angle types") != std::string::npos) ||
           (line.find("dihedral types") != std::string::npos) ||
           (line.find("improper types") != std::string::npos) ||
           (line.find("extra bond per atom") != std::string::npos) ||
           (line.find("extra angle per atom") != std::string::npos) ||
           (line.find("extra dihedral per atom") != std::string::npos) ||
           (line.find("extra improper per atom") != std::string::npos) ||
           (line.find("extra special per atom") != std::string::npos) ||
           (line.find("ellipsoids") != std::string::npos) ||
           (line.find("lines") != std::string::npos) ||
           (line.find("triangles") != std::string::npos) ||
           (line.find("bodies") != std::string::npos);
}

std::vector<size_t> guess_molecules(const Frame& frame) {
    // Initialize the molids vector with each atom in its own molecule
    auto molids = std::vector<size_t>();
    molids.reserve(frame.size());
    for (size_t i=0; i<frame.size(); i++) {
        molids.push_back(i);
    }

    for (auto bond: frame.topology().bonds()) {
        auto i = bond[0];
        auto j = bond[1];

        // Merge the bigger molid in the smaller one
        size_t new_id = molids[i];
        size_t old_id = molids[j];
        if (molids[i] > molids[j]) {
            new_id = molids[j];
            old_id = molids[i];
        }

        for (auto& molid: molids) {
            if (molid == old_id) {
                molid = new_id;
            }
        }
    }

    // Make sure the molids are consecutive
    std::unordered_map<size_t, size_t> molids_mapping;
    for (auto& molid: molids) {
        auto it = molids_mapping.find(molid);
        if (it != molids_mapping.end()) {
            molid = it->second;
        } else {
            // We've not found this id yet
            auto new_id = molids_mapping.size();
            molids_mapping.insert({molid, new_id});
            molid = new_id;
        }
    }

    return molids;
}

double tilt_factor(const Matrix3D& matrix, size_t i, size_t j) {
    assert(i != j);
    auto factor = matrix[i][j];
    auto a = matrix[i][i];
    if (factor >= 0) {
        while (fabs(factor) > a / 2) {
            factor -= a;
        }
    } else {
        while (factor < -a / 2) {
            factor += a;
        }
    }

    // Correct for rounding errors in calculations of sin/cos
    if (fabs(factor) < 1e-15) {
        factor = 0;
    }

    return factor;
}

optional<uint64_t> LAMMPSDataFormat::forward() {
    // LAMMPS Data only supports one step, so always act like there is only one
    auto position = file_.tellpos();
    if (position == 0) {
        // advance the pointer for the next call
        file_.readline();
        return position;
    } else {
        return nullopt;
    }
}
