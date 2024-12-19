// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cstddef>
#include <cstdint>
#include <vector>

#include "chemfiles/capi/types.h"
#include "chemfiles/capi/misc.h"
#include "chemfiles/capi/utils.hpp"
#include "chemfiles/capi/shared_allocator.hpp"

#include "chemfiles/capi/topology.h"

#include "chemfiles/Frame.hpp"
#include "chemfiles/Topology.hpp"
#include "chemfiles/Connectivity.hpp"

using namespace chemfiles;

extern "C" CHFL_TOPOLOGY* chfl_topology(void) {
    CHFL_TOPOLOGY* topology = nullptr;
    CHFL_ERROR_GOTO(
        topology = shared_allocator::make_shared<Topology>();
    )
    return topology;
error:
    chfl_free(topology);
    return nullptr;
}

extern "C" const CHFL_TOPOLOGY* chfl_topology_from_frame(const CHFL_FRAME* const frame) {
    const CHFL_TOPOLOGY* topology = nullptr;
    CHECK_POINTER_GOTO(frame);
    CHFL_ERROR_GOTO(
        topology = shared_allocator::shared_ptr<Topology>(frame, &frame->topology());
    )
    return topology;
error:
    chfl_free(topology);
    return nullptr;
}

extern "C" CHFL_TOPOLOGY* chfl_topology_copy(const CHFL_TOPOLOGY* const topology) {
    CHFL_TOPOLOGY* new_topology = nullptr;
    CHFL_ERROR_GOTO(
        new_topology = shared_allocator::make_shared<Topology>(*topology);
    )
    return new_topology;
error:
    chfl_free(new_topology);
    return nullptr;
}

extern "C" chfl_status chfl_topology_atoms_count(const CHFL_TOPOLOGY* const topology, uint64_t* const count) {
    CHECK_POINTER(topology);
    CHECK_POINTER(count);
    CHFL_ERROR_CATCH(
        *count = static_cast<uint64_t>(topology->size());
    )
}

extern "C" chfl_status chfl_topology_resize(CHFL_TOPOLOGY* const topology, uint64_t natoms) {
    CHECK_POINTER(topology);
    CHFL_ERROR_CATCH(
        topology->resize(checked_cast(natoms));
    )
}

extern "C" chfl_status chfl_topology_add_atom(CHFL_TOPOLOGY* const topology, const CHFL_ATOM* const atom) {
    CHECK_POINTER(topology);
    CHECK_POINTER(atom);
    CHFL_ERROR_CATCH(
        topology->add_atom(*atom);
    )
}

extern "C" chfl_status chfl_topology_remove(CHFL_TOPOLOGY* const topology, uint64_t i) {
    CHECK_POINTER(topology);
    CHFL_ERROR_CATCH(
        topology->remove(checked_cast(i));
    )
}

extern "C" chfl_status chfl_topology_bonds_count(const CHFL_TOPOLOGY* const topology, uint64_t* const count) {
    CHECK_POINTER(topology);
    CHECK_POINTER(count);
    CHFL_ERROR_CATCH(
        *count = static_cast<uint64_t>(topology->bonds().size());
    )
}

extern "C" chfl_status chfl_topology_angles_count(const CHFL_TOPOLOGY* const topology, uint64_t* const count) {
    CHECK_POINTER(topology);
    CHECK_POINTER(count);
    CHFL_ERROR_CATCH(
        *count = static_cast<uint64_t>(topology->angles().size());
    )
}

extern "C" chfl_status chfl_topology_dihedrals_count(const CHFL_TOPOLOGY* const topology, uint64_t* const count) {
    CHECK_POINTER(topology);
    CHECK_POINTER(count);
    CHFL_ERROR_CATCH(
        *count = static_cast<uint64_t>(topology->dihedrals().size());
    )
}

extern "C" chfl_status chfl_topology_impropers_count(const CHFL_TOPOLOGY* const topology, uint64_t* const count) {
    CHECK_POINTER(topology);
    CHECK_POINTER(count);
    CHFL_ERROR_CATCH(
        *count = static_cast<uint64_t>(topology->impropers().size());
    )
}

extern "C" chfl_status chfl_topology_bonds(const CHFL_TOPOLOGY* const topology, uint64_t (*data)[2], uint64_t count) {
    CHECK_POINTER(topology);
    CHECK_POINTER(data);

    CHFL_ERROR_CATCH(
        auto& bonds = topology->bonds();
        if (checked_cast(count) != bonds.size()) {
            set_last_error("wrong data size in function 'chfl_topology_bonds'.");
            return CHFL_MEMORY_ERROR;
        }

        for (size_t i=0; i<bonds.size(); i++) {
            data[i][0] = static_cast<uint64_t>(bonds[i][0]);
            data[i][1] = static_cast<uint64_t>(bonds[i][1]);
        }
    )
}

extern "C" chfl_status chfl_topology_angles(const CHFL_TOPOLOGY* const topology, uint64_t (*data)[3], uint64_t count) {
    CHECK_POINTER(topology);
    CHECK_POINTER(data);

    CHFL_ERROR_CATCH(
        auto& angles = topology->angles();
        if (checked_cast(count) != angles.size()) {
            set_last_error("wrong data size in function 'chfl_topology_angles'.");
            return CHFL_MEMORY_ERROR;
        }

        for (size_t i=0; i<angles.size(); i++) {
            data[i][0] = static_cast<uint64_t>(angles[i][0]);
            data[i][1] = static_cast<uint64_t>(angles[i][1]);
            data[i][2] = static_cast<uint64_t>(angles[i][2]);
        }
    )
}

extern "C" chfl_status chfl_topology_dihedrals(const CHFL_TOPOLOGY* const topology, uint64_t (*data)[4], uint64_t count) {
    CHECK_POINTER(topology);
    CHECK_POINTER(data);

    CHFL_ERROR_CATCH(
        auto& dihedrals = topology->dihedrals();
        if (checked_cast(count) != dihedrals.size()) {
            set_last_error("wrong data size in function 'chfl_topology_dihedrals'.");
            return CHFL_MEMORY_ERROR;
        }

        for (size_t i=0; i<dihedrals.size(); i++) {
            data[i][0] = static_cast<uint64_t>(dihedrals[i][0]);
            data[i][1] = static_cast<uint64_t>(dihedrals[i][1]);
            data[i][2] = static_cast<uint64_t>(dihedrals[i][2]);
            data[i][3] = static_cast<uint64_t>(dihedrals[i][3]);
        }
    )
}

extern "C" chfl_status chfl_topology_impropers(const CHFL_TOPOLOGY* const topology, uint64_t (*data)[4], uint64_t count) {
    CHECK_POINTER(topology);
    CHECK_POINTER(data);

    CHFL_ERROR_CATCH(
        auto& impropers = topology->impropers();
        if (checked_cast(count) != impropers.size()) {
            set_last_error("wrong data size in function 'chfl_topology_impropers'.");
            return CHFL_MEMORY_ERROR;
        }

        for (size_t i=0; i<impropers.size(); i++) {
            data[i][0] = static_cast<uint64_t>(impropers[i][0]);
            data[i][1] = static_cast<uint64_t>(impropers[i][1]);
            data[i][2] = static_cast<uint64_t>(impropers[i][2]);
            data[i][3] = static_cast<uint64_t>(impropers[i][3]);
        }
    )
}

extern "C" chfl_status chfl_topology_add_bond(CHFL_TOPOLOGY* const topology, uint64_t i, uint64_t j) {
    CHECK_POINTER(topology);
    CHFL_ERROR_CATCH(
        topology->add_bond(checked_cast(i), checked_cast(j));
    )
}

extern "C" chfl_status chfl_topology_remove_bond(CHFL_TOPOLOGY* const topology, uint64_t i, uint64_t j) {
    CHECK_POINTER(topology);
    CHFL_ERROR_CATCH(
        topology->remove_bond(checked_cast(i), checked_cast(j));
    )
}

extern "C" chfl_status chfl_topology_clear_bonds(CHFL_TOPOLOGY* const topology) {
    CHECK_POINTER(topology);
    CHFL_ERROR_CATCH(
        topology->clear_bonds();
    )
}

extern "C" chfl_status chfl_topology_residues_count(const CHFL_TOPOLOGY* const topology, uint64_t* const count) {
    CHECK_POINTER(topology);
    CHECK_POINTER(count);
    CHFL_ERROR_CATCH(
        *count = static_cast<uint64_t>(topology->residues().size());
    )
}

extern "C" chfl_status chfl_topology_add_residue(CHFL_TOPOLOGY* const topology, const CHFL_RESIDUE* const residue) {
    CHECK_POINTER(topology);
    CHECK_POINTER(residue);
    CHFL_ERROR_CATCH(
        topology->add_residue(*residue);
    )
}

extern "C" chfl_status chfl_topology_residues_linked(
    const CHFL_TOPOLOGY* const topology, const CHFL_RESIDUE* const first,
    const CHFL_RESIDUE* const second, bool* result) {
    CHECK_POINTER(topology);
    CHECK_POINTER(first);
    CHECK_POINTER(second);
    CHECK_POINTER(result);
    CHFL_ERROR_CATCH(
        *result = topology->are_linked(*first, *second);
    )
}

extern "C" chfl_status chfl_topology_bond_with_order(CHFL_TOPOLOGY* const topology, uint64_t i, uint64_t j, chfl_bond_order bond_order) {
    CHECK_POINTER(topology);
    CHFL_ERROR_CATCH(
        topology->add_bond(checked_cast(i), checked_cast(j), static_cast<Bond::BondOrder>(bond_order));
    )
}

extern "C" chfl_status chfl_topology_bond_orders(const CHFL_TOPOLOGY* const topology, chfl_bond_order orders[], uint64_t nbonds) {
    CHECK_POINTER(topology);
    CHECK_POINTER(orders);

    CHFL_ERROR_CATCH(
        if (nbonds != topology->bond_orders().size()) {
            set_last_error("wrong data size in function 'chfl_topology_bond_orders'.");
            return CHFL_MEMORY_ERROR;
        }

        auto& bond_orders = topology->bond_orders();
        for (size_t i=0; i<nbonds; i++) {
            orders[i] = static_cast<chfl_bond_order>(bond_orders[i]);
        }
    )
}

extern "C" chfl_status chfl_topology_bond_order( const CHFL_TOPOLOGY* const topology, uint64_t i, uint64_t j, chfl_bond_order* order) {
    CHECK_POINTER(topology);
    CHECK_POINTER(order);

    CHFL_ERROR_CATCH(
        *order = static_cast<chfl_bond_order>(
            topology->bond_order(checked_cast(i), checked_cast(j))
        );
    )
}
