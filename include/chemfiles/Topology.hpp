// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_TOPOLOGY_HPP
#define CHEMFILES_TOPOLOGY_HPP

#include <cstddef>
#include <string>
#include <vector>
#include <unordered_map>

#include "chemfiles/Atom.hpp"
#include "chemfiles/Connectivity.hpp"
#include "chemfiles/Residue.hpp"
#include "chemfiles/Error.hpp"
#include "chemfiles/exports.h"

#include "chemfiles/external/optional.hpp"

namespace chemfiles {

/// A topology contains the definition of all the atoms in the system, as well
/// as the liaisons between the particles (bonds, angles, dihedrals, ...) and
/// the residues.
///
/// Only the atoms and the bonds are stored, the angles, dihedrals and impropers
/// are automaticaly deduced from the bonds.
///
/// It is also possible to iterate over a `Topology`, yielding all the atoms in
/// the system.
///
/// @example{topology/iterate.cpp}
class CHFL_EXPORT Topology final {
public:
    using iterator = std::vector<Atom>::iterator;
    using const_iterator = std::vector<Atom>::const_iterator;

    /// Construct a new empty topology
    ///
    /// @example{topology/topology.cpp}
    Topology() = default;

    ~Topology() = default;
    Topology(const Topology&) = default;
    Topology& operator=(const Topology&) = default;
    Topology(Topology&&) = default;
    Topology& operator=(Topology&&) = default;

    /// Get a reference to the atom at the position `index`.
    ///
    /// @example{topology/index.cpp}
    ///
    /// @param index the atomic index
    /// @throws OutOfBounds if `index` is greater than `size()`
    Atom& operator[](size_t index) {
        if (index >= size()) {
            throw OutOfBounds(
                "atomic index out of bounds in topology: we have "
                + std::to_string(size()) + " atoms, but the index is "
                + std::to_string(index)
            );
        }
        return atoms_[index];
    }

    /// Get a const reference to the atom at the position `index`.
    ///
    /// @example{topology/index.cpp}
    ///
    /// @param index the atomic index
    /// @throws OutOfBounds if `index` is greater than `size()`
    const Atom& operator[](size_t index) const {
        if (index >= size()) {
            throw OutOfBounds(
                "atomic index out of bounds in topology: we have "
                + std::to_string(size()) + " atoms, but the index is "
                + std::to_string(index)
            );
        }
        return atoms_[index];
    }

    iterator begin() {return atoms_.begin();}
    const_iterator begin() const {return atoms_.begin();}
    const_iterator cbegin() const {return atoms_.cbegin();}
    iterator end() {return atoms_.end();}
    const_iterator end() const {return atoms_.end();}
    const_iterator cend() const {return atoms_.cend();}

    /// Add an `atom` at the end of this topology.
    ///
    /// @example{topology/add_atom.cpp}
    ///
    /// @param atom the new atom to add
    void add_atom(Atom atom);

    /// Delete the atom at index `i` in this topology, as well as all the bonds
    /// involving this atom.
    ///
    /// This function modify the index of all the atoms after `i`, and modify
    /// the bond list accordingly.
    ///
    /// @example{topology/remove.cpp}
    ///
    /// @param i the index of the atom to remove
    /// @throws OutOfBounds if `i` is greater than size()
    void remove(size_t i);

    /// Add a bond in the system, between the atoms at index `atom_i` and
    /// `atom_j`.
    ///
    /// @example{topology/add_bond.cpp}
    ///
    /// @param atom_i the index of the first atom in the bond
    /// @param atom_j the index of the second atom in the bond
    /// @param bond_order the bond order for the bond added
    /// @throws OutOfBounds if `atom_i` or `atom_j` are greater than `size()`
    /// @throws Error if `atom_i == atom_j`, as this is an invalid bond
    void add_bond(size_t atom_i, size_t atom_j, Bond::BondOrder bond_order = Bond::UNKNOWN);

    /// Remove a bond in the system, between the atoms at index `atom_i` and
    /// `atom_j`.
    ///
    /// If the bond does not exist, this does nothing.
    ///
    /// @example{topology/remove_bond.cpp}
    ///
    /// @param atom_i the index of the first atom in the bond
    /// @param atom_j the index of the second atom in the bond
    /// @throws OutOfBounds if `atom_i` or `atom_j` are greater than `size()`
    void remove_bond(size_t atom_i, size_t atom_j);

    /// Get the bond order for the given bond
    ///
    /// If the bond does not exist, this will thrown an Error.
    ///
    /// @example{topology/bond_order.cpp}
    ///
    /// @param atom_i the index of the first atom in the bond
    /// @param atom_j the index of the second atom in the bond
    /// @throws OutOfBounds if `atom_i` or `atom_j` are greater than `size()`
    /// @throws Error if no bond between `atom_i` and `atom_j` exists.
    Bond::BondOrder bond_order(size_t atom_i, size_t atom_j) const;

    /// Get the number of atoms in the topology
    ///
    /// @example{topology/size.cpp}
    size_t size() const {
        return atoms_.size();
    }

    /// Resize the topology to hold `size` atoms, adding new atoms as needed.
    ///
    /// If the new number of atoms is bigger than the old one, pre-existing
    /// atoms are conserved.
    ///
    /// If the new size if smaller than the old one, all atoms and connectivity
    /// elements after the new size are removed.
    ///
    /// @example{topology/resize.cpp}
    ///
    /// @param size the new size of the topology
    void resize(size_t size);

    /// Allocate memory in the frame to be able to store data for `size` atoms.
    ///
    /// This function does not change the actual number of atoms in the
    /// topology, and should be used as an optimisation.
    ///
    /// @example{topology/reserve.cpp}
    ///
    /// @param size the number of elements to reserve memory for
    void reserve(size_t size);

    /// Get the bonds in the system
    ///
    /// The bonds are sorted according to `operator<(const Bond&, const Bond&)`,
    /// which mean it is possible to look for a bond in the list using a binary
    /// search (`std::lower_bound`).
    ///
    /// @example{topology/bonds.cpp}
    const std::vector<Bond>& bonds() const;

    /// Get the bond orders in the system.
    ///
    /// The bond orders are sorted so that the index of each bond is the same as
    /// its index in the array returned by `Topology::bonds`. This means that
    /// the bond order for `Topology::bonds()[index]` would be given by
    /// `bond_orders()[index]`.
    ///
    /// @example{topology/bond_order.cpp}
    const std::vector<Bond::BondOrder>& bond_orders() const;

    /// Get the angles in the system
    ///
    /// The angles are sorted according to `operator<(const Angle&, const
    /// Angle&)`, which mean it is possible to look for an angle in the list
    /// using a binary search (`std::lower_bound`).
    ///
    /// @example{topology/angles.cpp}
    const std::vector<Angle>& angles() const;

    /// Get the dihedral angles in the system
    ///
    /// The dihedrals are sorted according to `operator<(const Dihedral&, const
    /// Dihedral&)`, which mean it is possible to look for a dihedral in the
    /// list using a binary search (`std::lower_bound`).
    ///
    /// @example{topology/dihedrals.cpp}
    const std::vector<Dihedral>& dihedrals() const;

    /// Get the improper dihedral angles in the system
    ///
    /// The impropers are sorted according to `operator<(const Improper&, const
    /// Improper&)`, which mean it is possible to look for an improper in the
    /// list using a binary search (`std::lower_bound`).
    ///
    /// @example{topology/impropers.cpp}
    const std::vector<Improper>& impropers() const;

    /// Remove all bonding information in the topology (bonds, angles and
    /// dihedrals)
    ///
    /// @example{topology/clear_bonds.cpp}
    void clear_bonds() {
        connect_ = Connectivity();
    }

    /// Add a `residue` to this topology.
    ///
    /// @example{topology/add_residue.cpp}
    ///
    /// @param residue the residue to add to this topology
    /// @throw chemfiles::Error if any atom in the `residue` is already in
    ///        another residue in this topology. In that case, the topology is
    ///        not modified.
    void add_residue(Residue residue);

    /// Check if two residues are linked together, i.e. if there is a bond
    /// between one atom in the `first` residue and one atom in the `second`
    /// one. Both residues should be in this topology.
    ///
    /// The two residues are the same (`first == second`), this function returns
    /// `true`.
    ///
    /// @example{topology/are_linked.cpp}
    bool are_linked(const Residue& first, const Residue& second) const;

    /// Get the residue containing the atom at the given `index`.
    ///
    /// If no residue contains this atom, this function returns `nullopt`.
    ///
    /// @example{topology/residue_for_atom.cpp}
    optional<const Residue&> residue_for_atom(size_t index) const;

    /// Get the residue at the given `index` in this topology
    ///
    /// There is no guarantee that this index matches the residue id.
    ///
    /// @example{topology/residue.cpp}
    const Residue& residue(size_t index) const {
        if (index >= residues_.size()) {
            throw OutOfBounds(
                "residue index out of bounds in topology: we have "
                + std::to_string(residues_.size()) + " residues, "
                + "but the index is " + std::to_string(index)
            );
        }
        return residues_[index];
    }

    /// Get all the residues in the topology as a vector
    ///
    /// @example{topology/residues.cpp}
    const std::vector<Residue>& residues() const {
        return residues_;
    }

private:
    /// Atoms in the system.
    std::vector<Atom> atoms_;
    /// Connectivity of the system.
    Connectivity connect_;
    /// List of residues in the system.
    std::vector<Residue> residues_;
    /// Association between atom indexes and residues indexes.
    std::unordered_map<size_t, size_t> residue_mapping_;
};

} // namespace chemfiles

#endif
