// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cstdint>
#include <cstring>
#include <string>

#include "chemfiles/capi/types.h"
#include "chemfiles/capi/misc.h"
#include "chemfiles/capi/utils.hpp"
#include "chemfiles/capi/shared_allocator.hpp"

#include "chemfiles/capi/residue.h"

#include "chemfiles/error_fmt.hpp"
#include "chemfiles/Property.hpp"
#include "chemfiles/Residue.hpp"
#include "chemfiles/Topology.hpp"  // IWYU pragma: keep

#include "chemfiles/external/optional.hpp"

using namespace chemfiles;

extern "C" CHFL_RESIDUE* chfl_residue(const char* name) {
    CHFL_RESIDUE* residue = nullptr;
    CHECK_POINTER_GOTO(name);
    CHFL_ERROR_GOTO(
        residue = shared_allocator::make_shared<Residue>(std::string(name));
    )
    return residue;
error:
    chfl_free(residue);
    return nullptr;
}

extern "C" CHFL_RESIDUE* chfl_residue_with_id(const char* name, int64_t resid) {
    CHFL_RESIDUE* residue = nullptr;
    CHECK_POINTER_GOTO(name);
    CHFL_ERROR_GOTO(
        residue = shared_allocator::make_shared<Residue>(std::string(name), resid);
    )
    return residue;
error:
    chfl_free(residue);
    return nullptr;
}

extern "C" const CHFL_RESIDUE* chfl_residue_from_topology(const CHFL_TOPOLOGY* const topology, uint64_t i) {
    const CHFL_RESIDUE* residue = nullptr;
    CHECK_POINTER_GOTO(topology);
    CHFL_ERROR_GOTO(
        residue = shared_allocator::shared_ptr<Residue>(topology, &topology->residue(checked_cast(i)));
    )
    return residue;
error:
    chfl_free(residue);
    return nullptr;
}

extern "C" const CHFL_RESIDUE* chfl_residue_for_atom(const CHFL_TOPOLOGY* const topology, uint64_t i) {
    const CHFL_RESIDUE* residue = nullptr;
    CHECK_POINTER_GOTO(topology);
    CHFL_ERROR_GOTO(
        auto optional = topology->residue_for_atom(checked_cast(i));
        if (optional) {
            residue = shared_allocator::shared_ptr<Residue>(topology, &*optional);
        }
    )
    return residue;
error:
    chfl_free(residue);
    return nullptr;
}

extern "C" CHFL_RESIDUE* chfl_residue_copy(const CHFL_RESIDUE* const residue) {
    CHFL_RESIDUE* new_residue = nullptr;
    CHFL_ERROR_GOTO(
        new_residue = shared_allocator::make_shared<Residue>(*residue);
    )
    return new_residue;
error:
    chfl_free(new_residue);
    return nullptr;
}

extern "C" chfl_status chfl_residue_atoms_count(const CHFL_RESIDUE* const residue, uint64_t* const count) {
    CHECK_POINTER(residue);
    CHECK_POINTER(count);
    CHFL_ERROR_CATCH(
        *count = static_cast<uint64_t>(residue->size());
    )
}

extern "C" chfl_status chfl_residue_atoms(const CHFL_RESIDUE* const residue, uint64_t atoms[], uint64_t count) {
    CHECK_POINTER(residue);
    CHECK_POINTER(atoms);
    CHFL_ERROR_CATCH(
        if (checked_cast(count) != residue->size()) {
            set_last_error("wrong data size in function 'chfl_residue_atoms'.");
            return CHFL_MEMORY_ERROR;
        }

        size_t i = 0;
        for (size_t index: *residue) {
            atoms[i] = static_cast<uint64_t>(index);
            i++;
        }
    )
}

extern "C" chfl_status chfl_residue_id(const CHFL_RESIDUE* const residue, int64_t* id) {
    CHECK_POINTER(residue);
    CHECK_POINTER(id);
    CHFL_ERROR_CATCH(
        try {
            *id = residue->id().value();
        } catch (const bad_optional_access&) {
            throw Error("this residue does not have a id");
        }
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

extern "C" chfl_status chfl_residue_properties_count(const CHFL_RESIDUE* const residue, uint64_t* const count) {
    CHECK_POINTER(residue);
    CHECK_POINTER(count);
    CHFL_ERROR_CATCH(
        *count = static_cast<uint64_t>(residue->properties().size());
    )
}

extern "C" chfl_status chfl_residue_list_properties(const CHFL_RESIDUE* const residue, const char* names[], uint64_t count) {
    CHECK_POINTER(residue);
    CHECK_POINTER(names);
    CHFL_ERROR_CATCH(
        auto& properties = residue->properties();
        if (checked_cast(count) != properties.size()) {
            set_last_error("wrong data size in function 'chfl_frame_list_properties'.");
            return CHFL_MEMORY_ERROR;
        }

        size_t i = 0;
        for (auto& it: properties) {
            names[i] = it.first.c_str();
            i++;
        }
    )
}

extern "C" chfl_status chfl_residue_set_property(CHFL_RESIDUE* const residue, const char* name, const CHFL_PROPERTY* const property) {
    CHECK_POINTER(residue);
    CHECK_POINTER(name);
    CHECK_POINTER(property);
    CHFL_ERROR_CATCH(
        residue->set(name, *property);
    )
}

extern "C" CHFL_PROPERTY* chfl_residue_get_property(const CHFL_RESIDUE* const residue, const char* name) {
    CHFL_PROPERTY* property = nullptr;
    CHECK_POINTER_GOTO(residue);
    CHECK_POINTER_GOTO(name);
    CHFL_ERROR_GOTO(
        auto residue_property = residue->get(name);
        if (residue_property) {
            property = shared_allocator::make_shared<Property>(*residue_property);
        } else {
            throw property_error("can not find a property named '{}' in this residue", name);
        }
    )
    return property;
error:
    chfl_free(property);
    return nullptr;
}
