// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>

#include "chemfiles/Residue.hpp"
#include "chemfiles/error_fmt.hpp"
#include "chemfiles/external/optional.hpp"

using namespace chemfiles;

Residue::Residue(std::string name): name_(std::move(name)), id_(nullopt) {}

Residue::Residue(std::string name, int64_t resid): name_(std::move(name)), id_(resid) {}

void Residue::add_atom(size_t i) {
    atoms_.insert(i);
}

bool Residue::contains(size_t i) const {
    return atoms_.find(i) != atoms_.end();
}

void Residue::remove(size_t i) {
    auto iter = atoms_.find(i);
    if (iter == atoms_.end()) {
        // we should only remove an atom if we know it exists in the residue
        throw error("invalid call to Residue::remove, this is a bug");
    }
    atoms_.erase(iter);
}

void Residue::atom_removed(size_t i) {
    for (auto& atom: atoms_.as_mutable_vec()) {
        if (atom > i) {
            atom -= 1;
        }
    }
}
