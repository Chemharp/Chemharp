// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_FRAME_HPP
#define CHEMFILES_FRAME_HPP

#include <cstddef>
#include <string>
#include <vector>

#include "chemfiles/exports.h"
#include "chemfiles/types.hpp"
#include "chemfiles/external/span.hpp"
#include "chemfiles/external/optional.hpp"

#include "chemfiles/Topology.hpp"
#include "chemfiles/UnitCell.hpp"
#include "chemfiles/Connectivity.hpp"
#include "chemfiles/Property.hpp"
#include "chemfiles/Residue.hpp"

namespace chemfiles {
class Atom;

/// A frame contains data from one simulation step The Frame class holds data
/// from one step of a simulation: the current topology, the positions, and the
/// velocities of the particles in the system. If some information is missing
/// the corresponding data is filled with a default value. Specifically:
///
/// @verbatim embed:rst:leading-slashes
///
/// * ``cell`` is an infinite unit cell;
/// * ``topology`` is empty, and contains no data;
/// * ``positions`` is filled with zeros;
/// * ``velocities`` is the ``nullopt`` variant of :cpp:class:`chemfiles::optional`.
///
/// @endverbatim
///
///Iterating over a `Frame` will yield all the atoms in the system.
///
/// @example{frame/iterate.cpp}
class CHFL_EXPORT Frame final {
public:
    /// Create an empty frame with no atoms and the given cell.
    ///
    /// @example{frame/frame.cpp}
    explicit Frame(UnitCell cell = UnitCell());

    ~Frame() = default;
    Frame(Frame&&) = default;
    Frame& operator=(Frame&&) = default;

    /// Get a clone (exact copy) of this frame.
    ///
    /// This replace the implicit copy constructor (which is private) to
    /// make an explicit copy of the frame.
    ///
    /// @example{frame/clone.cpp}
    Frame clone() const {
        return *this;
    }

    /// Get a const reference to the topology of this frame
    ///
    /// It is not possible to get a modifiable reference to the topology,
    /// because it would then be possible to remove/add atoms without changing
    /// the actual positions and velocity storage. Instead, all the mutating
    /// functionalities of the topology are mirrored on the frame (adding and
    /// removing bonds, adding residues, *etc.*)
    ///
    /// @example{frame/topology.cpp}
    const Topology& topology() const {
        return topology_;
    }

    /// Set the topology of this frame to `topology`
    ///
    /// @param topology the new Topology to use for this frame
    ///
    /// @throw Error if the topology size does not match the size of this frame
    ///
    /// @example{frame/topology.cpp}
    void set_topology(Topology topology);

    /// Get a const reference to the unit cell of this frame
    ///
    /// @example{frame/cell.cpp}
    const UnitCell& cell() const {
        return cell_;
    }

    /// Get a reference to the unit cell of this frame
    ///
    /// @example{frame/cell.cpp}
    UnitCell& cell() {
        return cell_;
    }

    /// Set the unit cell for this frame to `cell`
    ///
    /// @param cell the new UnitCell to use for this frame
    ///
    /// @example{frame/cell.cpp}
    void set_cell(UnitCell cell) {
        cell_ = std::move(cell);  // NOLINT: std::move for trivially copyable type
    }

    /// Get the number of atoms in this frame
    ///
    /// @example{frame/size.cpp}
    size_t size() const;

    /// Get the positions (in Angstroms) of the atoms in this frame.
    ///
    /// @example{frame/positions.cpp}
    span<Vector3D> positions() {
        return positions_;
    }

    /// Get the positions (in Angstroms) of the atoms in this frame as a const
    /// reference
    ///
    /// @example{frame/positions.cpp}
    const std::vector<Vector3D>& positions() const {
        return positions_;
    }

    /// Add velocities data storage to this frame.
    ///
    /// If velocities are already defined, this functions does nothing. The new
    /// velocities are initialized to 0.
    ///
    /// @example{frame/add_velocities.cpp}
    void add_velocities();

    /// Get an velocities (in Angstroms/ps) of the atoms in this frame, if this
    /// frame contains velocity data.
    ///
    /// @example{frame/velocities.cpp}
    optional<span<Vector3D>> velocities() {
        if (velocities_) {
            return {*velocities_};
        } else {
            return nullopt;
        }
    }

    /// Get an velocities (in Angstroms/ps) of the atoms in this frame as a
    /// const reference, if this frame contains velocity data.
    ///
    /// @example{frame/velocities.cpp}
    optional<const std::vector<Vector3D>&> velocities() const {
        if (velocities_) {
            return {*velocities_};
        } else {
            return nullopt;
        }
    }

    /// Resize the frame to contain `size` atoms.
    ///
    /// If the new number of atoms is bigger than the old one, missing data is
    /// initialized to 0. Pre-existing values are conserved.
    ///
    /// If the new size if smaller than the old one, all atoms and connectivity
    /// elements after the new size are removed.
    ///
    /// @example{frame/resize.cpp}
    void resize(size_t size);

    /// Allocate memory in the frame to have enough size for `size` atoms.
    ///
    /// This function does not change the actual number of atoms in the frame,
    /// and should be used as an optimisation.
    ///
    /// @example{frame/reserve.cpp}
    void reserve(size_t size);

    /// Add an `atom` at the given `position` and optionally with the given
    /// `velocity`. The `velocity` value will only be used if this frame
    /// contains velocity data.
    ///
    /// @example{frame/add_atom.cpp}
    void add_atom(Atom atom, Vector3D position, Vector3D velocity = Vector3D());

    /// Remove the atom at index `i` in the system.
    ///
    /// @throws chemfiles::OutOfBounds if `i` is bigger than the number of atoms
    ///         in this frame
    ///
    /// @example{frame/remove.cpp}
    void remove(size_t i);

    /// Get the current simulation step.
    ///
    /// The step is set by the `Trajectory` when reading a frame.
    ///
    /// @example{frame/step.cpp}
    size_t step() const {
        return step_;
    }

    /// Set the current simulation step to `step`
    ///
    /// @example{frame/step.cpp}
    void set_step(size_t step) {
        step_ = step;
    }

    /// Guess the bonds, angles, dihedrals and impropers angles in this frame.
    ///
    /// The bonds are guessed using a distance-based algorithm, and then angles,
    /// dihedrals and impropers are guessed from the bonds. The distance
    /// criterion uses the Van der Waals radii of the atoms.
    ///
    /// @throw Error if the Van der Waals radius in unknown for a given atom.
    ///
    /// @example{frame/guess_bonds.cpp}
    void guess_bonds();

    /// Remove all connectivity information in the frame's topology
    ///
    /// @example{frame/clear_bonds.cpp}
    void clear_bonds() {
        topology_.clear_bonds();
    }

    /// Add a `residue` to this frame's topology.
    ///
    /// @example{frame/add_residue.cpp}
    ///
    /// @param residue the residue to add to this topology
    /// @throw chemfiles::Error if any atom in the `residue` is already in
    ///        another residue in this topology. In that case, the topology is
    ///        not modified.
    void add_residue(Residue residue) {
        topology_.add_residue(std::move(residue));
    }

    /// Add a bond in the system, between the atoms at index `atom_i` and
    /// `atom_j`.
    ///
    /// @example{frame/add_bond.cpp}
    ///
    /// @param atom_i the index of the first atom in the bond
    /// @param atom_j the index of the second atom in the bond
    /// @param bond_order the bond order of the new bond
    /// @throws OutOfBounds if `atom_i` or `atom_j` are greater than `size()`
    /// @throws Error if `atom_i == atom_j`, as this is an invalid bond
    void add_bond(size_t atom_i, size_t atom_j, Bond::BondOrder bond_order = Bond::UNKNOWN, std::string bond_type = "") {
        topology_.add_bond(atom_i, atom_j, bond_order, bond_type);
    }

    /// Add an angle in the system, between the atoms at index `atom_i`,
    /// `atom_j` and `atom_k`.
    ///
    /// @example{frame/add_angle.cpp}
    ///
    /// @param atom_i the index of the first atom in the angle
    /// @param atom_j the index of the second atom in the angle
    /// @param atom_k the index of the third atom in the angle
    /// @param angle_type the angle order of the new angle
    /// @throws OutOfBounds if `atom_i`, `atom_j` and `atom_k` are 
    /// greater than `size()`
    /// @throws Error if any two indices are equal as this is an invalid angle
    void add_angle(size_t atom_i, size_t atom_j, size_t atom_k, std::string angle_type = "") {
        topology_.add_angle(atom_i, atom_j, atom_k, angle_type);
    }

    /// Add a dihedral in the system, between the atoms at index `atom_i`,
    /// `atom_j`, `atom_k`, and `atom_l`.
    ///
    /// @example{frame/add_dihedral.cpp}
    ///
    /// @param atom_i the index of the first atom in the dihedral
    /// @param atom_j the index of the second atom in the dihedral
    /// @param atom_k the index of the third atom in the dihedral
    /// @param atom_l the index of the forth atom in the dihedral
    /// @param dihedral_type the dihedral order of the new dihedral
    /// @throws OutOfBounds if `atom_i`, `atom_j` and `atom_k` are 
    /// greater than `size()`
    /// @throws Error if any two indices are equal as this is an invalid dihedral
    void add_dihedral(size_t atom_i, size_t atom_j, size_t atom_k, size_t atom_l, std::string dihedral_type = "") {
        topology_.add_dihedral(atom_i, atom_j, atom_k, atom_l, dihedral_type);
    }

    /// Add a improper in the system, between the atoms at index `atom_i`,
    /// `atom_j` and `atom_k`.
    ///
    /// @example{frame/add_improper.cpp}
    ///
    /// @param atom_i the index of the first atom in the improper
    /// @param atom_j the index of the second atom in the improper
    /// @param atom_k the index of the third atom in the improper
    /// @param atom_l the index of the forth atom in the improper
    /// @param improper_type the improper order of the new improper
    /// @throws OutOfBounds if `atom_i`, `atom_j` and `atom_k` are 
    /// greater than `size()`
    /// @throws Error if any two indices are equal as this is an invalid improper
    void add_improper(size_t atom_i, size_t atom_j, size_t atom_k, size_t atom_l, std::string improper_type = "") {
        topology_.add_improper(atom_i, atom_j, atom_k, atom_l, improper_type);
    }

    /// Remove a bond in the system, between the atoms at index `atom_i` and
    /// `atom_j`.
    ///
    /// If the bond does not exist, this does nothing.
    ///
    /// @example{frame/remove_bond.cpp}
    ///
    /// @param atom_i the index of the first atom in the bond
    /// @param atom_j the index of the second atom in the bond
    /// @throws OutOfBounds if `atom_i` or `atom_j` are greater than `size()`
    void remove_bond(size_t atom_i, size_t atom_j) {
        topology_.remove_bond(atom_i, atom_j);
    }

    /// Get a reference to the atom at the position `index`.
    ///
    /// @example{frame/index.cpp}
    ///
    /// @param index the atomic index
    /// @throws OutOfBounds if `index` is greater than `size()`
    Atom& operator[](size_t index) {
        return topology_[index];
    }

    /// Get a const reference to the atom at the position `index`.
    ///
    /// @example{frame/index.cpp}
    ///
    /// @param index the atomic index
    /// @throws OutOfBounds if `index` is greater than `size()`
    const Atom& operator[](size_t index) const {
        return topology_[index];
    }

    using iterator = Topology::iterator;
    using const_iterator = Topology::const_iterator;
    iterator begin() {return topology_.begin();}
    const_iterator begin() const {return topology_.begin();}
    const_iterator cbegin() const {return topology_.cbegin();}
    iterator end() {return topology_.end();}
    const_iterator end() const {return topology_.end();}
    const_iterator cend() const {return topology_.cend();}

    /// Get the distance between the atoms at indexes `i` and `j`, accounting
    /// for periodic boundary conditions. The distance is expressed in angstroms.
    ///
    /// @throws chemfiles::OutOfBounds if `i` or `j` are bigger than the number
    ///         of atoms in this frame
    ///
    /// @example{frame/distance.cpp}
    double distance(size_t i, size_t j) const;

    /// Get the angle formed by the atoms at indexes `i`, `j` and `k`,
    /// accounting for periodic boundary conditions. The angle is expressed in
    /// radians.
    ///
    /// @throws chemfiles::OutOfBounds if `i`, `j` or `k` are bigger than the
    ///         number of atoms in this frame
    ///
    /// @example{frame/angle.cpp}
    double angle(size_t i, size_t j, size_t k) const;

    /// Get the dihedral angle formed by the atoms at indexes `i`, `j`, `k` and
    /// `m`, accounting for periodic boundary conditions. The angle is expressed
    /// in radians.
    ///
    /// @throws chemfiles::OutOfBounds if `i`, `j`, `k` or `m` are bigger than
    ///         the number of atoms in this frame
    ///
    /// @example{frame/dihedral.cpp}
    double dihedral(size_t i, size_t j, size_t k, size_t m) const;

    /// Get the out of plane distance formed by the atoms at indexes `i`, `j`,
    /// `k` and `m`, accounting for periodic boundary conditions. The distance
    /// is expressed in angstroms.
    ///
    /// This is the distance between the atom j and the ikm plane. The j atom
    /// is the center of the improper dihedral angle formed by i, j, k and m.
    ///
    /// @throws chemfiles::OutOfBounds if `i`, `j`, `k` or `m` are bigger than
    ///         the number of atoms in this frame
    ///
    /// @example{frame/out_of_plane.cpp}
    double out_of_plane(size_t i, size_t j, size_t k, size_t m) const;

    /// Get the map of properties associated with this frame. This map might be
    /// iterated over to list the properties of the frame, or directly accessed.
    ///
    /// @example{frame/properties.cpp}
    const property_map& properties() const {
        return properties_;
    }

    /// Set an arbitrary property for this frame with the given `name` and
    /// `value`. If a property with this name already exist, it is silently
    /// replaced with the new value.
    ///
    /// @example{frame/set.cpp}
    void set(std::string name, Property value) {
        properties_.set(std::move(name), std::move(value));
    }

    /// Get the `Property` with the given `name` for this frame if it exists.
    ///
    /// If no property with the given `name` is found, this function returns
    /// `nullopt`.
    ///
    /// @example{frame/get.cpp}
    optional<const Property&> get(const std::string& name) const {
        return properties_.get(name);
    }

    /// Get the `Property` with the given `name` for this frame if it exists,
    /// and check that it has the required `kind`.
    ///
    /// If no property with the given `name` is found, this function returns
    /// `nullopt`.
    ///
    /// If a property with the given `name` is found, but has a different kind,
    /// this function emits a warning and returns `nullopt`.
    ///
    /// @example{frame/get.cpp}
    template<Property::Kind kind>
    optional<typename property_metadata<kind>::type> get(const std::string& name) const {
        return properties_.get<kind>(name);
    }

private:
    Frame(const Frame&) = default;
    Frame& operator=(const Frame&) = default;

    /// Current simulation step
    size_t step_ = 0;
    /// Positions of the particles
    std::vector<Vector3D> positions_;
    /// Velocities of the particles
    optional<std::vector<Vector3D>> velocities_;
    /// Topology of the described system
    Topology topology_;
    /// Unit cell of the system
    UnitCell cell_;
    /// Properties stored in this frame
    property_map properties_;
};

} // namespace chemfiles

#endif
