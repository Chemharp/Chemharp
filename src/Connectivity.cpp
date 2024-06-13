// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cassert>
#include <array>
#include <cstddef>
#include <vector>
#include <iterator>
#include <algorithm>

#include "chemfiles/Connectivity.hpp"
#include "chemfiles/error_fmt.hpp"
#include "chemfiles/sorted_set.hpp"

using namespace chemfiles;

Bond::Bond(size_t i, size_t j): data_({{0}}) {
    if (i == j) {
        throw error("can not have a bond between an atom and itself");
    }

    data_[0] = std::min(i, j);
    data_[1] = std::max(i, j);
}

size_t Bond::operator[](size_t i) const {
    if (i >= 2) {
        throw out_of_bounds("can not access atom n° {} in bond", i);
    }
    return data_[i];
}

Angle::Angle(size_t i, size_t j, size_t k): data_({{0}}) {
    if (i == j || i == k || j == k) {
        throw error("can not have the same atom twice in an angle");
    }

    data_[0] = std::min(i, k);
    data_[1] = j;
    data_[2] = std::max(i, k);
}

size_t Angle::operator[](size_t i) const {
    if (i >= 3) {
        throw out_of_bounds("can not access atom n° {} in angle", i);
    }
    return data_[i];
}

Dihedral::Dihedral(size_t i, size_t j, size_t k, size_t m): data_({{0}}) {
    if (i == j || j == k || k == m) {
        throw error("can not have an atom linked to itself in a dihedral angle");
    }

    if (i == k || j == m || i == m) {
        throw error("can not have an atom twice in a dihedral angle");
    }

    if (std::max(i, j) < std::max(k, m)) {
        data_[0] = i;
        data_[1] = j;
        data_[2] = k;
        data_[3] = m;
    } else {
        data_[0] = m;
        data_[1] = k;
        data_[2] = j;
        data_[3] = i;
    }
}

size_t Dihedral::operator[](size_t i) const {
    if (i >= 4) {
        throw out_of_bounds("can not access atom n° {} in dihedral", i);
    }
    return data_[i];
}

Improper::Improper(size_t i, size_t j, size_t k, size_t m): data_({{0}}) {
    if (j == i || j == k || j == m) {
        throw error("can not have an atom linked to itself in an improper dihedral angle");
    }

    if (i == k || i == m || k == m) {
        throw error("can not have an atom twice in an improper dihedral angle");
    }

    std::array<size_t, 3> others = {{i, k, m}};
    std::sort(others.begin(), others.end());
    data_[0] = others[0];
    data_[1] = j;
    data_[2] = others[1];
    data_[3] = others[2];
}

size_t Improper::operator[](size_t i) const {
    if (i >= 4) {
        throw out_of_bounds("can not access atom n° {} in improper", i);
    }
    return data_[i];
}

void Connectivity::recalculate() const {
    angles_.clear();
    dihedrals_.clear();
    impropers_.clear();
    angle_types_.clear();
    dihedral_types_.clear();
    improper_types_.clear();

    // Generate the list of which atom is bonded to which one
    auto bonded_to = std::vector<std::vector<size_t>>(biggest_atom_ + 1);
    for (auto const& bond: bonds_) {
        assert(bond[0] < bonded_to.size());
        assert(bond[1] < bonded_to.size());
        bonded_to[bond[0]].push_back(bond[1]);
        bonded_to[bond[1]].push_back(bond[0]);
    }

    // Generate the list of angles
    for (auto const& bond: bonds_) {
        auto i = bond[0];
        auto j = bond[1];
        for (auto k: bonded_to[i]) {
            if (k != j) {
                angles_.insert(Angle(k, i, j));
                angle_types_.push_back("");
            }
        }

        for (auto k: bonded_to[j]) {
            if (k != i) {
                angles_.insert(Angle(i, j, k));
                angle_types_.push_back("");
            }
        }
    }

    // Generate the list of dihedrals
    for (auto const& angle: angles_) {
        auto i = angle[0];
        auto j = angle[1];
        auto k = angle[2];
        for (auto m: bonded_to[i]) {
            if (m != j && m != k) {
                dihedrals_.insert(Dihedral(m, i, j, k));
                dihedral_types_.push_back("");
            }
        }

        for (auto m: bonded_to[k]) {
            if (m != i && m != j) {
                dihedrals_.insert(Dihedral(i, j, k, m));
                dihedral_types_.push_back("");
            }
        }

        for (auto m: bonded_to[j]) {
            if (m != i && m != k) {
                impropers_.insert(Improper(i, j, k, m));
                improper_types_.push_back("");
            }
        }
    }
    uptodate_ = true;
}

const sorted_set<Bond>& Connectivity::bonds() const {
    return bonds_;
}

const std::vector<Bond::BondOrder>& Connectivity::bond_orders() const {
    return bond_orders_;
}

const sorted_set<Angle>& Connectivity::angles() const {
    if (!uptodate_) {
        recalculate();
    }
    return angles_;
}

const sorted_set<Dihedral>& Connectivity::dihedrals() const {
    if (!uptodate_) {
        recalculate();
    }
    return dihedrals_;
}

const sorted_set<Improper>& Connectivity::impropers() const {
    if (!uptodate_) {
        recalculate();
    }
    return impropers_;
}

void Connectivity::add_bond(size_t i, size_t j, Bond::BondOrder bond_order, std::string bond_type) {
    uptodate_ = false;
    auto result = bonds_.emplace(i, j);
    if (i > biggest_atom_) {biggest_atom_ = i;}
    if (j > biggest_atom_) {biggest_atom_ = j;}

    if (result.second) {
        auto diff = std::distance(bonds_.cbegin(), result.first);
        bond_orders_.insert(bond_orders_.begin() + diff, bond_order);
        bond_types_.insert(bond_types_.begin() + diff, std::move(bond_type));
    }
}

void Connectivity::add_angle(size_t i, size_t j, size_t k, std::string angle_type) {
    uptodate_ = true;
    auto result = angles_.emplace(i, j, k);
    if (i > biggest_atom_) {biggest_atom_ = i;}
    if (j > biggest_atom_) {biggest_atom_ = j;}
    if (k > biggest_atom_) {biggest_atom_ = k;}

    if (result.second) {
        auto diff = std::distance(angles_.cbegin(), result.first);
        angle_types_.insert(angle_types_.begin() + diff, std::move(angle_type));
    }
}

void Connectivity::add_dihedral(size_t i, size_t j, size_t k, size_t l, std::string dihedral_type) {
    uptodate_ = true;
    auto result = dihedrals_.emplace(i, j, k, l);
    if (i > biggest_atom_) {biggest_atom_ = i;}
    if (j > biggest_atom_) {biggest_atom_ = j;}
    if (k > biggest_atom_) {biggest_atom_ = k;}
    if (l > biggest_atom_) {biggest_atom_ = l;}

    if (result.second) {
        auto diff = std::distance(dihedrals_.cbegin(), result.first);
        dihedral_types_.insert(dihedral_types_.begin() + diff, std::move(dihedral_type));
    }
}

void Connectivity::add_improper(size_t i, size_t j, size_t k, size_t l, std::string improper_type) {
    uptodate_ = true;
    auto result = impropers_.emplace(i, j, k, l);
    if (i > biggest_atom_) {biggest_atom_ = i;}
    if (j > biggest_atom_) {biggest_atom_ = j;}
    if (k > biggest_atom_) {biggest_atom_ = k;}
    if (l > biggest_atom_) {biggest_atom_ = l;}

    if (result.second) {
        auto diff = std::distance(impropers_.cbegin(), result.first);
        improper_types_.insert(improper_types_.begin() + diff, std::move(improper_type));
    }
}

void Connectivity::remove_bond(size_t i, size_t j) {
    auto pos = bonds_.find(Bond(i, j));
    if (pos != bonds_.end()) {
        uptodate_ = false;
        auto result = bonds_.erase(pos);

        auto diff = std::distance(bonds_.cbegin(), result);
        bond_orders_.erase(bond_orders_.begin() + diff);
        assert(bond_orders_.size() == bonds_.size());
    }
}

void Connectivity::atom_removed(size_t index) {
    auto to_remove = std::vector<Bond>();
    auto to_add = std::vector<Bond>();
    auto bo_add = std::vector<Bond::BondOrder>();

    const auto& bonds = bonds_.as_vec();
    for (size_t idx = 0; idx < bonds_.size(); idx++) {
        const auto& bond = bonds[idx];

        if (bond[0] == index || bond[1] == index) {
            throw error("can not shift atomic indexes that still have a bond");
        }
        if (bond[0] > index || bond[1] > index) {
            to_remove.push_back(bond);

            auto i = bond[0] > index ? bond[0] - 1 : bond[0];
            auto j = bond[1] > index ? bond[1] - 1 : bond[1];
            to_add.emplace_back(i, j);
            bo_add.push_back(bond_orders_[idx]);
        }
    }

    for (auto bond: to_remove) {
        this->remove_bond(bond[0], bond[1]);
    }

    for (size_t idx = 0; idx < to_add.size(); idx++) {
        const auto& bond = to_add[idx];
        this->add_bond(bond[0], bond[1], bo_add[idx]);
    }
}

Bond::BondOrder Connectivity::bond_order(size_t i, size_t j) const {
    auto pos = bonds_.find(Bond(i, j));
    if (pos != bonds_.end()) {
        auto diff = std::distance(bonds_.cbegin(), pos);
        return bond_orders_[static_cast<size_t>(diff)];
    }

    throw error(
        "out of bounds atomic index in `Connectivity::bond_order`: "
        "No bond between {} and {} exists",
        i, j
    );
}

const std::string& Connectivity::bond_type(size_t i, size_t j) const {
    auto pos = bonds_.find(Bond(i, j));
    if (pos != bonds_.end()) {
        auto diff = std::distance(bonds_.cbegin(), pos);
        return bond_types_[static_cast<size_t>(diff)];
    }

    throw error(
        "out of bounds atomic index in `Connectivity::bond_order`: "
        "No bond between {} and {} exists",
        i, j
    );
}

const std::string& Connectivity::angle_type(size_t i, size_t j, size_t k) const {
    auto pos = angles_.find(Angle(i, j, k));
    if (pos != angles_.end()) {
        auto diff = std::distance(angles_.cbegin(), pos);
        return angle_types_[static_cast<size_t>(diff)];
    }

    throw error(
        "out of bounds atomic index in `Connectivity::angle_type`: "
        "No angle between {}, {} and {} exists",
        i, j, k
    );
}

const std::string& Connectivity::dihedral_type(size_t i, size_t j, size_t k, size_t m) const {
    auto pos = dihedrals_.find(Dihedral(i, j, k, m));
    if (pos != dihedrals_.end()) {
        auto diff = std::distance(dihedrals_.cbegin(), pos);
        return dihedral_types_[static_cast<size_t>(diff)];
    }

    throw error(
        "out of bounds atomic index in `Connectivity::dihedral_type`: "
        "No dihedral between {}, {}, {} and {} exists",
        i, j, k, m
    );
}

const std::string& Connectivity::improper_type(size_t i, size_t j, size_t k, size_t m) const {
    auto pos = impropers_.find(Improper(i, j, k, m));
    if (pos != impropers_.end()) {
        auto diff = std::distance(impropers_.cbegin(), pos);
        return improper_types_[static_cast<size_t>(diff)];
    }

    throw error(
        "out of bounds atomic index in `Connectivity::improper_type`: "
        "No improper between {}, {}, {} and {} exists",
        i, j, k, m
    );
}