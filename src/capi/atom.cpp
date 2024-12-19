// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cstdint>
#include <cstring>
#include <string>

#include "chemfiles/capi/types.h"
#include "chemfiles/capi/misc.h"
#include "chemfiles/capi/utils.hpp"
#include "chemfiles/capi/shared_allocator.hpp"

#include "chemfiles/capi/atom.h"

#include "chemfiles/Atom.hpp"
#include "chemfiles/Frame.hpp"
#include "chemfiles/error_fmt.hpp"
#include "chemfiles/Property.hpp"

#include "chemfiles/external/optional.hpp"

using namespace chemfiles;

extern "C" CHFL_ATOM* chfl_atom(const char* name) {
    CHFL_ATOM* atom = nullptr;
    CHFL_ERROR_GOTO(
        atom = shared_allocator::make_shared<Atom>(name);
    )
    return atom;
error:
    chfl_free(atom);
    return nullptr;
}

extern "C" CHFL_ATOM* chfl_atom_copy(const CHFL_ATOM* const atom) {
    CHFL_ATOM* new_atom = nullptr;
    CHFL_ERROR_GOTO(
        new_atom = shared_allocator::make_shared<Atom>(*atom);
    )
    return new_atom;
error:
    chfl_free(new_atom);
    return nullptr;
}

extern "C" CHFL_ATOM* chfl_atom_from_frame(CHFL_FRAME* const frame, uint64_t index) {
    CHFL_ATOM* atom = nullptr;
    CHECK_POINTER_GOTO(frame);
    CHFL_ERROR_GOTO(
        // Return NULL if the index is out of bounds
        if (index >= frame->size()) {
            throw out_of_bounds(
                "out of bounds atomic index in `chfl_atom_from_frame`: we have {} atoms, but the index is {}",
                frame->size(), index
            );
        }
        atom = shared_allocator::shared_ptr<Atom>(frame, &(*frame)[checked_cast(index)]);
    )
    return atom;
error:
    chfl_free(atom);
    return nullptr;
}

extern "C" CHFL_ATOM* chfl_atom_from_topology(CHFL_TOPOLOGY* const topology, uint64_t index) {
    CHFL_ATOM* atom = nullptr;
    CHECK_POINTER_GOTO(topology);
    CHFL_ERROR_GOTO(
        // Return NULL if the index is out of bounds
        if (index >= topology->size()) {
            throw out_of_bounds(
                "out of bounds atomic index in `chfl_atom_from_topology`: we have {} atoms, but the index is {}",
                topology->size(), index
            );
        }
        atom = shared_allocator::shared_ptr<Atom>(topology, &(*topology)[checked_cast(index)]);
    )
    return atom;
error:
    chfl_free(atom);
    return nullptr;
}

extern "C" chfl_status chfl_atom_mass(const CHFL_ATOM* const atom, double* mass) {
    CHECK_POINTER(atom);
    CHECK_POINTER(mass);
    CHFL_ERROR_CATCH(
        *mass = atom->mass();
    )
}

extern "C" chfl_status chfl_atom_set_mass(CHFL_ATOM* const atom, double mass) {
    CHECK_POINTER(atom);
    CHFL_ERROR_CATCH(
        atom->set_mass(mass);
    )
}

extern "C" chfl_status chfl_atom_charge(const CHFL_ATOM* const atom, double* charge) {
    CHECK_POINTER(atom);
    CHECK_POINTER(charge);
    CHFL_ERROR_CATCH(
        *charge = atom->charge();
    )
}

extern "C" chfl_status chfl_atom_set_charge(CHFL_ATOM* const atom, double charge) {
    CHECK_POINTER(atom);
    CHFL_ERROR_CATCH(
        atom->set_charge(charge);
    )
}

extern "C" chfl_status chfl_atom_type(const CHFL_ATOM* const atom, char* const type, uint64_t buffsize) {
    CHECK_POINTER(atom);
    CHECK_POINTER(type);
    CHFL_ERROR_CATCH(
        strncpy(type, atom->type().c_str(), checked_cast(buffsize) - 1);
        type[buffsize - 1] = '\0';
    )
}

extern "C" chfl_status chfl_atom_set_type(CHFL_ATOM* const atom, const char* type) {
    CHECK_POINTER(atom);
    CHECK_POINTER(type);
    CHFL_ERROR_CATCH(
        atom->set_type(type);
    )
}

extern "C" chfl_status chfl_atom_name(const CHFL_ATOM* const atom, char* const name, uint64_t buffsize) {
    CHECK_POINTER(atom);
    CHECK_POINTER(name);
    CHFL_ERROR_CATCH(
        strncpy(name, atom->name().c_str(), checked_cast(buffsize) - 1);
        name[buffsize - 1] = '\0';
    )
}

extern "C" chfl_status chfl_atom_set_name(CHFL_ATOM* const atom, const char* name) {
    CHECK_POINTER(atom);
    CHECK_POINTER(name);
    CHFL_ERROR_CATCH(
        atom->set_name(name);
    )
}

extern "C" chfl_status chfl_atom_full_name(const CHFL_ATOM* const atom, char* const name, uint64_t buffsize) {
    CHECK_POINTER(atom);
    CHECK_POINTER(name);
    CHFL_ERROR_CATCH(
        auto full_name = atom->full_name();
        if (full_name) {
            std::strncpy(name, full_name.value().c_str(), checked_cast(buffsize) - 1);
            name[buffsize - 1] = '\0';
        } else {
            std::memset(name, 0, checked_cast(buffsize));
        }
    )
}

extern "C" chfl_status chfl_atom_vdw_radius(const CHFL_ATOM* const atom, double* radius) {
    CHECK_POINTER(atom);
    CHECK_POINTER(radius);
    CHFL_ERROR_CATCH(
        *radius = atom->vdw_radius().value_or(0);
    )
}

extern "C" chfl_status chfl_atom_covalent_radius(const CHFL_ATOM* const atom, double* radius) {
    CHECK_POINTER(atom);
    CHECK_POINTER(radius);
    CHFL_ERROR_CATCH(
        *radius = atom->covalent_radius().value_or(0);
    )
}

extern "C" chfl_status chfl_atom_atomic_number(const CHFL_ATOM* const atom, uint64_t* number) {
    CHECK_POINTER(atom);
    CHECK_POINTER(number);
    CHFL_ERROR_CATCH(
        *number = atom->atomic_number().value_or(0ul);
    )
}

extern "C" chfl_status chfl_atom_properties_count(const CHFL_ATOM* const atom, uint64_t* const count) {
    CHECK_POINTER(atom);
    CHECK_POINTER(count);
    CHFL_ERROR_CATCH(
        if (atom->properties()) {
            *count = static_cast<uint64_t>(atom->properties().value().size());
        }
        else {
            *count = 0;
        }
    )
}

extern "C" chfl_status chfl_atom_list_properties(const CHFL_ATOM* const atom, const char* names[], uint64_t count) {
    CHECK_POINTER(atom);
    CHECK_POINTER(names);
    CHFL_ERROR_CATCH(
        auto& properties = atom->properties();
        size_t property_size = 0;
        if (properties) {
            property_size = (*properties).size();
        }
        if (checked_cast(count) != property_size) {
            set_last_error("wrong data size in function 'chfl_atom_list_properties'.");
            return CHFL_MEMORY_ERROR;
        }

        if (properties) {
            size_t i = 0;
            for (auto& it: *properties) {
                names[i] = it.first.c_str();
                i++;
            }
        }
    )
}

extern "C" chfl_status chfl_atom_set_property(CHFL_ATOM* const atom, const char* name, const CHFL_PROPERTY* const property) {
    CHECK_POINTER(atom);
    CHECK_POINTER(name);
    CHECK_POINTER(property);
    CHFL_ERROR_CATCH(
        atom->set(name, *property);
    )
}

extern "C" CHFL_PROPERTY* chfl_atom_get_property(const CHFL_ATOM* const atom, const char* name) {
    CHFL_PROPERTY* property = nullptr;
    CHECK_POINTER_GOTO(atom);
    CHECK_POINTER_GOTO(name);
    CHFL_ERROR_GOTO(
        auto atom_property = atom->get(name);
        if (atom_property) {
            property = shared_allocator::make_shared<Property>(*atom_property);
        } else {
            throw property_error("can not find a property named '{}' in this atom", name);
        }
    )
    return property;
error:
    chfl_free(property);
    return nullptr;
}
