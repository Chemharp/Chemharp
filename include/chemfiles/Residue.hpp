// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_RESIDUE_HPP
#define CHEMFILES_RESIDUE_HPP

#include <cstddef>
#include <cstdint>

#include <string>
#include <algorithm>
#include <utility>

#include "chemfiles/exports.h"
#include "chemfiles/sorted_set.hpp"
#include "chemfiles/external/optional.hpp"
#include "chemfiles/Property.hpp"

namespace chemfiles {

/// A `Residue` is a group of atoms belonging to the same logical unit. They can
/// be molecules, amino-acids in a protein, monomers in polymers, *etc.*
///
/// Iterating over a `Residue` gives the indexes of the atoms in the associated
/// `Topology`.
///
/// @example{residue/iterate.cpp}
class CHFL_EXPORT Residue final {
public:
    /// Create a new residue with a given `name` and no residue id.
    ///
    /// @example{residue/residue-1.cpp}
    explicit Residue(std::string name);

    /// Create a new residue with a given `name` and residue id `resid`.
    ///
    /// @example{residue/residue-2.cpp}
    Residue(std::string name, int64_t resid);

    ~Residue() = default;
    Residue(const Residue&) = default;
    Residue& operator=(const Residue&) = default;
    Residue(Residue&&) = default;
    Residue& operator=(Residue&&) = default;

    /// Get the name of the residue.
    ///
    /// @example{residue/name.cpp}
    const std::string& name() const {
        return name_;
    }

    /// Get the residue identifier, or `nullopt` if it does not exist.
    ///
    /// @example{residue/id.cpp}
    optional<int64_t> id() const {
        return id_;
    }

    /// Get the size of the residue, i.e. the number of atoms in this residue.
    ///
    /// @example{residue/size.cpp}
    size_t size() const {
        return atoms_.size();
    }

    /// Add an atom with index `i` to this residue
    ///
    /// If the atom is already in the residue, this does nothing.
    ///
    /// @example{residue/add_atom.cpp}
    void add_atom(size_t i);

    /// Check if the residue contains a given atom with index `i`
    ///
    /// @example{residue/contains.cpp}
    bool contains(size_t i) const;

    using const_iterator = sorted_set<size_t>::const_iterator;
    // Iterators over the indexes of the atoms in the residue
    const_iterator begin() const {return atoms_.begin();}
    const_iterator end() const {return atoms_.end();}
    const_iterator cbegin() const {return atoms_.cbegin();}
    const_iterator cend() const {return atoms_.cend();}

    /// Get the map of properties asociated with this residue. This map might be
    /// iterated over to list the properties of the residue, or directly
    /// accessed.
    ///
    /// @example{frame/properties.cpp}
    const property_map& properties() const {
        return properties_;
    }

    /// Set an arbitrary `Property` for this residue with the given `name` and
    /// `value`. If a property with this name already exist, it is replaced with
    /// the new value.
    ///
    /// @example{residue/property.cpp}
    void set(std::string name, Property value) {
        properties_.set(std::move(name), std::move(value));
    }

    /// Get the `Property` with the given `name` for this residue if it exists.
    ///
    /// If no property with the given `name` is found, this function returns
    /// `nullopt`.
    ///
    /// @example{residue/property.cpp}
    optional<const Property&> get(const std::string& name) const {
        return properties_.get(name);
    }

    /// Get the `Property` with the given `name` for this residue if it exists,
    /// and check that it has the required `kind`.
    ///
    /// If no property with the given `name` is found, this function returns
    /// `nullopt`.
    ///
    /// If a property with the given `name` is found, but has a different kind,
    /// this function emits a warning and returns `nullopt`.
    ///
    /// @example{residue/property.cpp}
    template<Property::Kind kind>
    optional<typename property_metadata<kind>::type> get(const std::string& name) const {
        return properties_.get<kind>(name);
    }

private:
    /// Name of the residue
    std::string name_;
    /// Index of the residue in the initial topology file
    optional<int64_t> id_;
    /// Indexes of the atoms in this residue. These indexes refers to the
    /// associated topology.
    sorted_set<size_t> atoms_;
    /// Additional properties of this residue
    property_map properties_;

    /// Update the atomic indexes in this residue after an atom has been
    /// removed from the containing topology.
    ///
    /// This function shifts all the indexes bigger than `i` by -1.
    void atom_removed(size_t i);

    /// remove the atom at index i from this residue
    void remove(size_t i);

    friend bool operator==(const Residue& lhs, const Residue& rhs);

    friend class Topology;
};

inline bool operator==(const Residue& lhs, const Residue& rhs) {
    return lhs.id() == rhs.id() &&
           lhs.name() == rhs.name() &&
           lhs.size() == rhs.size() &&
           std::equal(lhs.begin(), lhs.end(), rhs.begin()) &&
           lhs.properties_ == rhs.properties_;
}

inline bool operator!=(const Residue& lhs, const Residue& rhs) {
    return !(lhs == rhs);
}

} // namespace chemfiles

#endif
