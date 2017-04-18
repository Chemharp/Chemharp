// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cstring>

#include "chemfiles/capi/residue.h"
#include "chemfiles/capi.hpp"

#include "chemfiles/Residue.hpp"
#include "chemfiles/Topology.hpp"
using namespace chemfiles;

extern "C" CHFL_RESIDUE* chfl_residue(const char* name, uint64_t resid) {
    CHFL_RESIDUE* residue = nullptr;
    CHECK_POINTER_GOTO(name);
    CHFL_ERROR_GOTO(
        residue = new Residue(std::string(name), resid);
    )
    return residue;
error:
    delete residue;
    return nullptr;
}

extern "C" CHFL_RESIDUE* chfl_residue_from_topology(const CHFL_TOPOLOGY* const topology, uint64_t i) {
    CHFL_RESIDUE* residue = nullptr;
    CHECK_POINTER_GOTO(topology);
    CHFL_ERROR_GOTO(
        // Return NULL if the index is out of bounds
        if (i >= topology->residues().size()) {
            throw MemoryError(
                "Out of bounds residue index " + std::to_string(i) +
                ". Last residue is " + std::to_string(topology->residues().size()) + "."
            );
        }
        residue = new Residue(topology->residues()[checked_cast(i)]);
    )
    return residue;
error:
    delete residue;
    return nullptr;
}

extern "C" CHFL_RESIDUE* chfl_residue_for_atom(const CHFL_TOPOLOGY* const topology, uint64_t i) {
    CHFL_RESIDUE* residue = nullptr;
    CHECK_POINTER_GOTO(topology);
    CHFL_ERROR_GOTO(
        auto res = topology->residue(checked_cast(i));
        if (res) {
            residue = new Residue(*res);
        }
    )
    return residue;
error:
    delete residue;
    return nullptr;
}

extern "C" CHFL_RESIDUE* chfl_residue_copy(const CHFL_RESIDUE* const residue) {
    CHFL_RESIDUE* new_residue = nullptr;
    CHFL_ERROR_GOTO(
        new_residue = new Residue(*residue);
    )
    return new_residue;
error:
    delete new_residue;
    return nullptr;
}

extern "C" chfl_status chfl_residue_atoms_count(const CHFL_RESIDUE* const residue, uint64_t* size) {
    CHECK_POINTER(residue);
    CHECK_POINTER(size);
    CHFL_ERROR_CATCH(
        *size = residue->size();
    )
}


extern "C" chfl_status chfl_residue_id(const CHFL_RESIDUE* const residue, uint64_t* id) {
    CHECK_POINTER(residue);
    CHECK_POINTER(id);
    CHFL_ERROR_CATCH(
        *id = residue->id();
    )
}

extern "C" chfl_status chfl_residue_name(const CHFL_RESIDUE* const residue, char* name, uint64_t buffsize) {
    CHECK_POINTER(residue);
    CHECK_POINTER(name);
    CHFL_ERROR_CATCH(
        strncpy(name, residue->name().c_str(), checked_cast(buffsize) - 1);
        name[buffsize - 1] = '\0';
    )
}

extern "C" chfl_status chfl_residue_add_atom(CHFL_RESIDUE* const residue, uint64_t i) {
    CHECK_POINTER(residue);
    CHFL_ERROR_CATCH(
        residue->add_atom(checked_cast(i));
    )
}

extern "C" chfl_status chfl_residue_contains(const CHFL_RESIDUE* const residue, uint64_t i, bool* result) {
    CHECK_POINTER(residue);
    CHECK_POINTER(result);
    CHFL_ERROR_CATCH(
        *result = residue->contains(checked_cast(i));
    )
}

extern "C" chfl_status chfl_residue_free(CHFL_RESIDUE* const residue) {
    delete residue;
    return CHFL_SUCCESS;
}
