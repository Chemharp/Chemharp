// // Chemfiles, a modern library for chemistry file reading and writing
// // Copyright (C) Guillaume Fraux and contributors -- BSD license

// #ifndef CHEMFILES_FORMAT_LAMMPSMol_HPP
// #define CHEMFILES_FORMAT_LAMMPSMol_HPP

// #include <cstdint>
// #include <string>
// #include <memory>

// #include "chemfiles/File.hpp"
// #include "chemfiles/Format.hpp"

// #include "chemfiles/external/optional.hpp"

// namespace chemfiles {
// class Frame;
// class MemoryBuffer;
// class FormatMetadata;

// class DataTypes {
// public:
//     DataTypes(const Topology& topology = Topology());

//     const sorted_set<atom_type>& atoms() const {return atoms_;}
//     const sorted_set<bond_type>& bonds() const {return bonds_;}
//     const sorted_set<angle_type>& angles() const {return angles_;}
//     const sorted_set<dihedral_type>& dihedrals() const {return dihedrals_;}
//     const sorted_set<improper_type>& impropers() const {return impropers_;}

//     /// Get the atom type number for the given atom.
//     ///
//     /// The atom must be in the topology used to construct this `DataTypes`
//     /// instance. The index numbering starts at zero, and can be used to index
//     /// the vector backing the `sorted_set<atom_type>` returned by `atoms()`.
//     size_t atom_type_id(const Atom& atom) const;

//     /// Get the bond type number for the bond type i-j.
//     ///
//     /// The bond type must be in the topology used to construct this `DataTypes`
//     /// instance. The index numbering starts at zero, and can be used to index
//     /// the vector backing the `sorted_set<bond_type>` returned by `bonds()`.
//     size_t bond_type_id(size_t type_i, size_t type_j) const;

//     /// Get the angle type number for the angle type i-j-k.
//     ///
//     /// The angle type must be in the topology used to construct this `DataTypes`
//     /// instance. The index numbering starts at zero, and can be used to index
//     /// the vector backing the `sorted_set<angle_type>` returned by `angles()`.
//     size_t angle_type_id(size_t type_i, size_t type_j, size_t type_k) const;

//     /// Get the dihedral type number for the dihedral type i-j-k-m.
//     ///
//     /// The dihedral type must be in the topology used to construct this
//     /// `DataTypes` instance. The index numbering starts at zero, and can be
//     /// used to index the vector backing the `sorted_set<dihedral_type>`
//     /// returned by `dihedrals()`.
//     size_t dihedral_type_id(size_t type_i, size_t type_j, size_t type_k, size_t type_m) const;

//     /// Get the improper type number for the improper type i-j-k-m.
//     ///
//     /// The improper type must be in the topology used to construct this
//     /// `DataTypes` instance. The index numbering starts at zero, and can be
//     /// used to index the vector backing the `sorted_set<improper_type>`
//     /// returned by `impropers()`.
//     size_t improper_type_id(size_t type_i, size_t type_j, size_t type_k, size_t type_m) const;

// private:
//     sorted_set<atom_type> atoms_;
//     sorted_set<bond_type> bonds_;
//     sorted_set<angle_type> angles_;
//     sorted_set<dihedral_type> dihedrals_;
//     sorted_set<improper_type> impropers_;
// };


// /// LAMMPS molecule template format reader and writer.
// ///
// class LAMMPSMolFormat final: public TextFormat {
// public:
//     LAMMPSMolFormat(std::string path, File::Mode mode, File::Compression compression):
//         TextFormat(std::move(path), mode, compression), current_section_(HEADER){}

//     LAMMPSMolFormat(std::shared_ptr<MemoryBuffer> memory, File::Mode mode, File::Compression compression) :
//         TextFormat(std::move(memory), mode, compression), current_section_(HEADER){}

//     void read_next(Frame& frame) override;
//     void write_next(const Frame& frame) override;
//     optional<uint64_t> forward() override;

// private:
//     enum section_t {
//         HEADER,
//         COORDS,
//         TYPES,
//         MOLECULES,
//         FRAGMENTS,
//         CHARGES,
//         DIAMETERS,
//         MASSES,
//         BONDS,
//         ANGLES,
//         DIHEDRALS,
//         IMPROPERS,
//         SPECIAL_BONDS_COUNTS,
//         SPECIAL_BONDS,
//         SHAKE_FLAGS,
//         SHAKE_ATOMS,
//         NOT_A_SECTION,
//         IGNORED;
//     } current_section_;

//     /// Get the section corresponding to a given line
//     section_t get_section(string_view line);

//     /// Read the header section
//     void read_header(Frame& frame);
//     size_t read_header_integer(string_view line, const std::string& context);

//     /// Get the section name from the next non-empty line
//     void get_next_section();
//     /// Skip all lines that are not sections names, and get the next section
//     void skip_to_next_section();

//     /// Read the atoms section
//     void read_coords(Frame& frame);
//     /// Read the masses section
//     void read_masses(Frame& frame);
//     /// Read the bonds section
//     void read_bonds(Frame& frame);

//     /// Write the header
//     void write_header(const Frame& frame);

//     size_t natoms_ = 0;
//     size_t nbonds_ = 0;

// };

// template<> const FormatMetadata& format_metadata<LAMMPSMolFormat>();

// } // namespace chemfiles

// #endif
