// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_FORMAT_LAMMPS_DATA_HPP
#define CHEMFILES_FORMAT_LAMMPS_DATA_HPP

#include <cstdint>
#include <tuple>
#include <limits>
#include <vector>
#include <string>
#include <memory>
#include <unordered_map>

#include "chemfiles/File.hpp"
#include "chemfiles/Format.hpp"

#include "chemfiles/Topology.hpp"  // IWYU pragma: keep
#include "chemfiles/sorted_set.hpp"
#include "chemfiles/external/optional.hpp"

namespace chemfiles {
class Atom;
class Frame;
class MemoryBuffer;
class FormatMetadata;

struct atom_data final {
    double x = 0;
    double y = 0;
    double z = 0;
    double charge = std::numeric_limits<double>::quiet_NaN();
    double mass = std::numeric_limits<double>::quiet_NaN();
    size_t index = 0;
    size_t type = 0;
    size_t molid = 0;
};

/// Possible LAMMPS atom style
class atom_style final {
private:
    /// The atom style name
    std::string name_;
    /// List of possible style for LAMMPS data files
    enum {
        ANGLE, ATOMIC, BODY, BOND, CHARGE, DIPOLE, DPD, ELECTRON, ELLIPSOID,
        FULL, LINE, MESO, MOLECULAR, PERI, SMD, SPHERE, TEMPLATE, TRI,
        WAVEPACKET, HYBRID
    } style_;
    /// Did we send the hybrid style warning?
    mutable bool warned_ = false;

public:
    explicit atom_style(std::string name);
    /// Read a single line with this atom style
    atom_data read_line(std::string_view line, size_t index) const;
};

// atom types are defined by the type string and the mass of the atom
using atom_type = std::pair<std::string, double>;
using bond_type = std::tuple<size_t, size_t>;
using angle_type = std::tuple<size_t, size_t, size_t>;
using dihedral_type = std::tuple<size_t, size_t, size_t, size_t>;
using improper_type = std::tuple<size_t, size_t, size_t, size_t>;

class DataTypes {
public:
    DataTypes(const Topology& topology = Topology());

    const sorted_set<atom_type>& atoms() const {return atoms_;}
    const sorted_set<std::tuple<bond_type, std::string>>& bonds() const {return bonds_;}
    const sorted_set<std::tuple<angle_type, std::string>>& angles() const {return angles_;}
    const sorted_set<std::tuple<dihedral_type, std::string>>& dihedrals() const {return dihedrals_;}
    const sorted_set<std::tuple<improper_type, std::string>>& impropers() const {return impropers_;}

    /// Get the atom type number for the given atom.
    ///
    /// The atom must be in the topology used to construct this `DataTypes`
    /// instance. The index numbering starts at zero, and can be used to index
    /// the vector backing the `sorted_set<atom_type>` returned by `atoms()`.
    size_t atom_type_id(const Atom& atom) const;

    /// Get the bond type number for the bond type i-j.
    ///
    /// The bond type must be in the topology used to construct this `DataTypes`
    /// instance. The index numbering starts at zero, and can be used to index
    /// the vector backing the `sorted_set<bond_type>` returned by `bonds()`.
    std::string bond_type_id(size_t type_i, size_t type_j) const;

    /// Get the angle type number for the angle type i-j-k.
    ///
    /// The angle type must be in the topology used to construct this `DataTypes`
    /// instance. The index numbering starts at zero, and can be used to index
    /// the vector backing the `sorted_set<angle_type, std::string>` returned by `angles()`.
    std::string angle_type_id(size_t type_i, size_t type_j, size_t type_k) const;

    /// Get the dihedral type number for the dihedral type i-j-k-m.
    ///
    /// The dihedral type must be in the topology used to construct this
    /// `DataTypes` instance. The index numbering starts at zero, and can be
    /// used to index the vector backing the `sorted_set<dihedral_type, std::string>`
    /// returned by `dihedrals()`.
    std::string dihedral_type_id(size_t type_i, size_t type_j, size_t type_k, size_t type_m) const;

    /// Get the improper type number for the improper type i-j-k-m.
    ///
    /// The improper type must be in the topology used to construct this
    /// `DataTypes` instance. The index numbering starts at zero, and can be
    /// used to index the vector backing the `sorted_set<improper_type, std::string>`
    /// returned by `impropers()`.
    std::string improper_type_id(size_t type_i, size_t type_j, size_t type_k, size_t type_m) const;

private:
    sorted_set<atom_type> atoms_;
    sorted_set<std::tuple<bond_type, std::string>> bonds_;
    sorted_set<std::tuple<angle_type, std::string>> angles_;
    sorted_set<std::tuple<dihedral_type, std::string>> dihedrals_;
    sorted_set<std::tuple<improper_type, std::string>> impropers_;
};

/// LAMMPS Data file format reader and writer.
///
/// LAMMPS data files are not fully stand-alone, as one needs to know the atom
/// style to read the data. This reader will try to guess the atom style by
/// checking the first line of the header for `atom_style <style>`, and by
/// reading any comment after the `Atoms` section name. If no atom style is
/// specified, the code default to `full` and send a warning.
///
/// The code also tries to read atomic names at the end of data lines. For
/// example, the atom at index 44 will have `C2` as atomic name.
///
/// ```
/// 44 44 2 0.000000 1.094000 2.061000 69.552002 # C2 RES
/// ```
class LAMMPSDataFormat final: public TextFormat {
public:
    LAMMPSDataFormat(std::string path, File::Mode mode, File::Compression compression):
        TextFormat(std::move(path), mode, compression),
        current_section_(HEADER),
        style_("full")
    {}

    LAMMPSDataFormat(std::shared_ptr<MemoryBuffer> memory, File::Mode mode, File::Compression compression):
        TextFormat(std::move(memory), mode, compression),
        current_section_(HEADER),
        style_("full")
    {}

    void read_next(Frame& frame) override;
    void write_next(const Frame& frame) override;
    optional<uint64_t> forward() override;

private:
    enum section_t {
        HEADER,
        ATOMS,
        MASSES,
        BONDS,
        ANGLES,
        DIHEDRALS,
        IMPROPERS,
        VELOCITIES,
        IGNORED,
        NOT_A_SECTION,
    } current_section_;

    /// Get the section corresponding to a given line
    section_t get_section(std::string_view line);

    /// Read the header section
    void read_header(Frame& frame);
    size_t read_header_integer(std::string_view line, const std::string& context);
    double read_header_box_bounds(std::string_view line, const std::string& lo, const std::string& hi);

    /// Get the section name from the next non-empty line
    void get_next_section();
    /// Skip all lines that are not sections names, and get the next section
    void skip_to_next_section();

    /// Read the atoms section
    void read_atoms(Frame& frame);
    /// Read the masses section
    void read_masses();
    /// Read the bonds section
    void read_bonds(Frame& frame);
    /// Read the angles section
    void read_angles(Frame& frame);
    /// Read the dihedrals section
    void read_dihedrals(Frame& frame);
    /// Read the impropers section
    void read_impropers(Frame& frame);
    /// Read the velocities section
    void read_velocities(Frame& frame);
    /// Setup masses of the frame with previously read values. This function
    /// must be called before `setup_names`.
    void setup_masses(Frame& frame) const;
    /// Setup atomic names with previously read values
    void setup_names(Frame& frame) const;

    /// Write the header
    void write_header(const DataTypes& types, const Frame& frame);
    /// Write the types sections
    void write_types(const DataTypes& types);
    /// Write the masses section
    void write_masses(const DataTypes& types);
    /// Write the Atoms section
    void write_atoms(const DataTypes& types, const Frame& frame);
    /// Write the Velocities section
    void write_velocities(const Frame& frame);
    /// Write the Bonds section
    void write_bonds(const DataTypes& types, const Topology& topology);
    /// Write the Angles section
    void write_angles(const DataTypes& types, const Topology& topology);
    /// Write the Dihedrals section
    void write_dihedrals(const DataTypes& types, const Topology& topology);
    /// Write the Impropers section
    void write_impropers(const DataTypes& types, const Topology& topology);

    // =============== Data used for reading files
    /// Name of the atom style to use when reading the "Atoms" section
    std::string atom_style_name_ = "";
    /// Atom style to use when reading the "Atoms" section
    atom_style style_;
    /// Number of atoms in the file
    size_t natoms_ = 0;
    /// Number of atom types in the file
    size_t natom_types_ = 0;
    /// Number of bonds in the file
    size_t nbonds_ = 0;
    /// Number of bond types in the file
    /// size_t nbond_types_ = 0;
    /// Number of angles in the file
    size_t nangles_ = 0;
    /// Number of angle types in the file
    /// size_t nangle_types_ = 0;
    /// Number of dihedrals in the file
    size_t ndihedrals_ = 0;
    /// Number of dihedral types in the file
    /// size_t ndihedral_types_ = 0;
    /// Number of impropers in the file
    size_t nimpropers_ = 0;
    /// Optional masses, indexed by atomic type
    std::unordered_map<std::string, double> masses_;
    /// Optional atomic names, indexed by atomic indexes
    std::vector<std::string> names_;
};

template<> const FormatMetadata& format_metadata<LAMMPSDataFormat>();

} // namespace chemfiles

#endif
