// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cstddef>
#include <utility>
#include <vector>
#include <unordered_map>

#include "chemfiles/Atom.hpp"
#include "chemfiles/Connectivity.hpp"
#include "chemfiles/error_fmt.hpp"
#include "chemfiles/Residue.hpp"
#include "chemfiles/Topology.hpp"
#include "chemfiles/sorted_set.hpp"
#include "chemfiles/external/optional.hpp"

using namespace chemfiles;

void Topology::resize(size_t size) {
    for (const auto& bond: connect_.bonds()) {
        if (bond[0] >= size || bond[1] >= size) {
            throw error(
                "can not resize the topology to contains {} atoms as there "
                "is a bond between atoms {} - {}",
                size, bond[0], bond[1]
            );
        }
    }
    atoms_.resize(size, Atom());
}

void Topology::add_atom(Atom atom) {
    atoms_.emplace_back(std::move(atom));
}

void Topology::reserve(size_t size) {
    atoms_.reserve(size);
}

void Topology::add_bond(size_t atom_i, size_t atom_j, Bond::BondOrder bond_order) {
    if (atom_i >= size() || atom_j >= size()) {
        throw out_of_bounds(
            "out of bounds atomic index in `Topology::add_bond`: "
            "we have {} atoms, but the bond indexes are {} and {}",
            size(), atom_i, atom_j
        );
    }
    connect_.add_bond(atom_i, atom_j, bond_order);
}

void Topology::remove_bond(size_t atom_i, size_t atom_j) {
    if (atom_i >= size() || atom_j >= size()) {
        throw out_of_bounds(
            "out of bounds atomic index in `Topology::remove_bond`: "
            "we have {} atoms, but the bond indexes are {} and {}",
            size(), atom_i, atom_j
        );
    }
    connect_.remove_bond(atom_i, atom_j);
}

Bond::BondOrder Topology::bond_order(size_t atom_i, size_t atom_j) const {
    if (atom_i >= size() || atom_j >= size()) {
        throw out_of_bounds(
            "out of bounds atomic index in `Topology::bond_order`: "
            "we have {} atoms, but the bond indexes are {} and {}",
            size(), atom_i, atom_j
        );
    }

    return connect_.bond_order(atom_i, atom_j);
}

void Topology::remove(size_t i) {
    if (i >= size()) {
        throw out_of_bounds(
            "out of bounds atomic index in `Topology::remove`: we have {} atoms, "
            "but the indexe is {}",
            size(), i
        );
    }
    atoms_.erase(atoms_.begin() + static_cast<std::ptrdiff_t>(i));

    // Remove all bonds with the removed atom
    auto bonds = connect_.bonds();
    for (const auto& bond : bonds) {
        if (bond[0] == i || bond[1] == i) {
            connect_.remove_bond(bond[0], bond[1]);
        }
    }
    // remove the atom from the corresponding residue
    auto it = residue_mapping_.find(i);
    if (it != residue_mapping_.end()) {
        residues_[it->second].remove(i);
    }

    // shift all bonds indexes
    connect_.atom_removed(i);
    // shift all residue atoms
    for (auto& res : residues_) {
        res.atom_removed(i);
    }
}

const std::vector<Bond>& Topology::bonds() const {
    return connect_.bonds().as_vec();
}

const std::vector<Bond::BondOrder>& Topology::bond_orders() const {
    return connect_.bond_orders();
}

const std::vector<Angle>& Topology::angles() const {
    return connect_.angles().as_vec();
}

const std::vector<Dihedral>& Topology::dihedrals() const {
    return connect_.dihedrals().as_vec();
}

const std::vector<Improper>& Topology::impropers() const {
    return connect_.impropers().as_vec();
}

void Topology::add_residue(Residue residue) {
    for (auto i: residue) {
        auto it = residue_mapping_.find(i);
        if (it != residue_mapping_.end()) {
            throw error(
                "can not add this residue: atom {} is already in another residue",
                i
            );
        }
    }
    auto res_index = residues_.size();
    residues_.emplace_back(std::move(residue));
    for (auto i: residues_.back()) {
        residue_mapping_.insert({i, res_index});
    }
}

bool Topology::are_linked(const Residue& first, const Residue& second) const {
    if (first == second) {
        return true;
    }
    auto bonds = connect_.bonds();
    for (auto i: first) {
        for (auto j: second) {
            if (bonds.find({i, j}) != bonds.end()) {
                return true;
            }
        }
    }
    return false;
}

optional<const Residue&> Topology::residue_for_atom(size_t index) const {
    auto it = residue_mapping_.find(index);
    if (it == residue_mapping_.end()) {
        // This atom is not in a residue
        return nullopt;
    } else {
        return residues_[it->second];
    }
}
