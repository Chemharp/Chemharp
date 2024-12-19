// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cassert>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <utility>

#include "chemfiles/Atom.hpp"
#include "chemfiles/periodic_table.hpp"

#include "chemfiles/utils.hpp"
#include "chemfiles/external/optional.hpp"

using namespace chemfiles;

static std::string normalize_atomic_name(const std::string& type) {
    assert(type.length() <= 2);
    std::string normalized = type;
    if (type.length() == 1) {
        normalized[0] = to_ascii_uppercase(normalized[0]);
    } else if (type.length() == 2) {
        normalized[0] = to_ascii_uppercase(normalized[0]);
        normalized[1] = to_ascii_lowercase(normalized[1]);
    }
    return normalized;
}

optional<const AtomicData&> chemfiles::find_in_periodic_table(const std::string& type) {
    atomic_data_map::const_iterator it;
    if (type.length() <= 2) {
        it = PERIODIC_TABLE.find(normalize_atomic_name(type));
    } else {
        it = PERIODIC_TABLE.find(type);
    }

    if (it != PERIODIC_TABLE.end()) {
        return it->second;
    } else {
        return nullopt;
    }
}

static optional<const AtomicData&> find_element(const std::string& type) {
    return find_in_periodic_table(type);
}

Atom::Atom(std::string name): name_(std::move(name)), type_(name_) {
    auto element = find_element(type_);
    if (element) {
        mass_ = element->mass.value_or(0);
        charge_ = element->charge.value_or(0);
    }
}

Atom::Atom(std::string name, std::string type): name_(std::move(name)), type_(std::move(type)) {
    auto element = find_element(type_);
    if (element) {
        mass_ = element->mass.value_or(0);
        charge_ = element->charge.value_or(0);
    }
}

optional<std::string> Atom::full_name() const {
    auto element = find_element(type_);
    if (element) {
        return element->full_name;
    } else {
        return nullopt;
    }
}

optional<double> Atom::vdw_radius() const {
    auto element = find_element(type_);
    if (element) {
        return element->vdw_radius;
    } else {
        return nullopt;
    }
}

optional<double> Atom::covalent_radius() const {
    auto element = find_element(type_);
    if (element) {
        return element->covalent_radius;
    } else {
        return nullopt;
    }
}

optional<uint64_t> Atom::atomic_number() const {
    auto element = find_element(type_);
    if (element) {
        return element->number;
    } else {
        return nullopt;
    }
}
