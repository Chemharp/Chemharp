// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cassert>
#include <cstddef>
#include <cmath>
#include <string>
#include <utility>
#include <vector>
#include <iterator>
#include <algorithm>
#include <unordered_map>

#include "chemfiles/types.hpp"
#include "chemfiles/error_fmt.hpp"
#include "chemfiles/external/optional.hpp"

#include "chemfiles/Atom.hpp"
#include "chemfiles/Topology.hpp"
#include "chemfiles/UnitCell.hpp"
#include "chemfiles/Connectivity.hpp"

#include "chemfiles/Frame.hpp"

using namespace chemfiles;

// get radius compatible with VMD bond guessing algorithm
static optional<double> guess_bonds_radius(const Atom& atom);

Frame::Frame(UnitCell cell): cell_(std::move(cell)) {} // NOLINT: std::move for trivially copyable type

size_t Frame::size() const {
    assert(positions_.size() == topology_.size());
    if (velocities_) {
        assert(positions_.size() == velocities_->size());
    }
    return positions_.size();
}

void Frame::resize(size_t size) {
    topology_.resize(size);
    positions_.resize(size);
    if (velocities_) {
        velocities_->resize(size);
    }
}

void Frame::reserve(size_t size) {
    topology_.reserve(size);
    positions_.reserve(size);
    if (velocities_) {
        velocities_->reserve(size);
    }
}

void Frame::add_velocities() {
    if (!velocities_) {
        velocities_ = std::vector<Vector3D>(size());
    }
}

void Frame::guess_bonds() {
    topology_.clear_bonds();
    // This bond guessing algorithm comes from VMD
    auto cutoff = 0.833;
    for (size_t i = 0; i < size(); i++) {
        auto rad = guess_bonds_radius(topology_[i]).value_or(0);
        cutoff = std::max(cutoff, rad);
    }
    cutoff = 1.2 * cutoff;

    for (size_t i = 0; i < size(); i++) {
        auto i_radius = guess_bonds_radius(topology_[i]);
        if (!i_radius) {
            throw error(
                "missing Van der Waals radius for '{}'", topology_[i].type()
            );
        }
        for (size_t j = i + 1; j < size(); j++) {
            auto j_radius = guess_bonds_radius(topology_[j]);
            if (!j_radius) {
                throw error(
                    "missing Van der Waals radius for '{}'", topology_[j].type()
                );
            }
            auto d = distance(i, j);
            auto radii = i_radius.value() + j_radius.value();
            if (0.03 < d && d < 0.6 * radii && d < cutoff) {
                topology_.add_bond(i, j);
            }
        }
    }

    auto bonds = topology().bonds();
    auto to_remove = std::vector<Bond>();
    // We need to remove bonds between hydrogen atoms which are bonded more than
    // once
    for (auto& bond : bonds) {
        auto i = bond[0];
        auto j = bond[1];
        if (topology_[i].type() != "H") {
            continue;
        }
        if (topology_[j].type() != "H") {
            continue;
        }

        auto nbonds = std::count_if(
            std::begin(bonds), std::end(bonds), [=](const Bond& b) {
                return b[0] == i || b[0] == j || b[1] == i || b[1] == j;
            });
        assert(nbonds >= 1);

        if (nbonds != 1) {
            to_remove.push_back(bond);
        }
    }

    for (auto& bond : to_remove) {
        topology_.remove_bond(bond[0], bond[1]);
    }
}

void Frame::set_topology(Topology topology) {
    if (topology.size() != size()) {
        throw error(
            "the topology contains {} atoms, but the frame contains {} atoms",
            topology.size(), size()
        );
    }
    topology_ = std::move(topology);
}

void Frame::add_atom(Atom atom, Vector3D position, Vector3D velocity) {
    topology_.add_atom(std::move(atom));
    positions_.push_back(position);
    if (velocities_) {
        velocities_->push_back(velocity);
    }
    assert(size() == topology_.size());
}

void Frame::remove(size_t i) {
    if (i >= size()) {
        throw out_of_bounds(
            "out of bounds atomic index in `Frame::remove`: we have {} atoms, "
            "but the index is {}",
            size(), i
        );
    }
    topology_.remove(i);
    positions_.erase(positions_.begin() + static_cast<std::ptrdiff_t>(i));
    if (velocities_) {
        velocities_->erase(velocities_->begin() + static_cast<std::ptrdiff_t>(i));
    }
    assert(size() == topology_.size());
}

double Frame::distance(size_t i, size_t j) const {
    if (i >= size() || j >= size()) {
        throw out_of_bounds(
            "out of bounds atomic index in `Frame::distance`: we have {} "
            "atoms, but the index are {} and {}",
            size(), i, j
        );
    }

    auto rij = positions_[i] - positions_[j];
    return cell_.wrap(rij).norm();
}

double Frame::angle(size_t i, size_t j, size_t k) const {
    if (i >= size() || j >= size() || k >= size()) {
        throw out_of_bounds(
            "out of bounds atomic index in `Frame::angle`: we have {} atoms, "
            "but the index are {}, {}, and {}",
            size(), i, j, k
        );
    }

    auto rij = cell_.wrap(positions_[i] - positions_[j]);
    auto rkj = cell_.wrap(positions_[k] - positions_[j]);

    auto cos = dot(rij, rkj) / (rij.norm() * rkj.norm());
    cos = std::max(-1.0, std::min(1.0, cos));
    return acos(cos);
}

double Frame::dihedral(size_t i, size_t j, size_t k, size_t m) const {
    if (i >= size() || j >= size() || k >= size() || m >= size()) {
        throw out_of_bounds(
            "out of bounds atomic index in `Frame::dihedral`: we have {} "
            "atoms, but the index are {}, {}, {}, and {}",
            size(), i, j, k, m
        );
    }

    auto rij = cell_.wrap(positions_[i] - positions_[j]);
    auto rjk = cell_.wrap(positions_[j] - positions_[k]);
    auto rkm = cell_.wrap(positions_[k] - positions_[m]);

    auto a = cross(rij, rjk);
    auto b = cross(rjk, rkm);
    return atan2(rjk.norm() * dot(b, rij), dot(a, b));
}

double Frame::out_of_plane(size_t i, size_t j, size_t k, size_t m) const {
    if (i >= size() || j >= size() || k >= size() || m >= size()) {
        throw out_of_bounds(
            "out of bounds atomic index in `Frame::out_of_plane`: we have {} "
            "atoms, but the index are {}, {}, {}, and {}",
            size(), i, j, k, m
        );
    }

    auto rji = cell_.wrap(positions_[j] - positions_[i]);
    auto rik = cell_.wrap(positions_[i] - positions_[k]);
    auto rim = cell_.wrap(positions_[i] - positions_[m]);

    auto n = cross(rik, rim);
    auto n_norm = n.norm();
    if (n_norm < 1e-12) {
        // if i, k, and m are colinear, then j is always inside the plane
        return 0;
    } else {
        return dot(rji, n) / n_norm;
    }
}

const std::unordered_map<std::string, double> BOND_GUESSING_RADII = {
    {"H", 1.0},
    {"C", 1.5},
    {"O", 1.3},
    {"N", 1.4},
    {"S", 1.9},
    {"F", 1.2},
};

optional<double> guess_bonds_radius(const Atom& atom) {
    const auto& type = atom.type();
    auto it = BOND_GUESSING_RADII.find(type);
    if (it != BOND_GUESSING_RADII.end()) {
        return it->second;
    } else {
        // default to chemfiles provided atom type
        return atom.vdw_radius();
    }
}
