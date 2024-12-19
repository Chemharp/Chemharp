// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cstddef>
#include <cstdint>
#include <string>

#include "chemfiles/capi/types.h"
#include "chemfiles/capi/misc.h"
#include "chemfiles/capi/utils.hpp"
#include "chemfiles/capi/shared_allocator.hpp"

#include "chemfiles/capi/frame.h"

#include "chemfiles/Frame.hpp"
#include "chemfiles/error_fmt.hpp"
#include "chemfiles/Property.hpp"
#include "chemfiles/Connectivity.hpp"

#include "chemfiles/types.hpp"
#include "chemfiles/external/optional.hpp"
#include "chemfiles/external/span.hpp"

using namespace chemfiles;

extern "C" CHFL_FRAME* chfl_frame(void) {
    CHFL_FRAME* frame = nullptr;
    CHFL_ERROR_GOTO(
        frame = shared_allocator::make_shared<Frame>();
    )
    return frame;
error:
    chfl_free(frame);
    return nullptr;
}

extern "C" CHFL_FRAME* chfl_frame_copy(const CHFL_FRAME* const frame) {
    CHFL_FRAME* new_frame = nullptr;
    CHFL_ERROR_GOTO(
        new_frame = shared_allocator::make_shared<Frame>(frame->clone());
    )
    return new_frame;
error:
    chfl_free(new_frame);
    return nullptr;
}

extern "C" chfl_status chfl_frame_atoms_count(const CHFL_FRAME* const frame, uint64_t* const count) {
    CHECK_POINTER(frame);
    CHECK_POINTER(count);
    CHFL_ERROR_CATCH(
        *count = static_cast<uint64_t>(frame->size());
    )
}

extern "C" chfl_status chfl_frame_positions(CHFL_FRAME* const frame, chfl_vector3d** positions, uint64_t* size) {
    CHECK_POINTER(frame);
    CHECK_POINTER(positions);
    CHECK_POINTER(size);
    static_assert(
        sizeof(chfl_vector3d) == sizeof(Vector3D),
        "Wrong size for chfl_vector3d. It should match Vector3D."
    );
    CHFL_ERROR_CATCH(
        auto cpp_positions = frame->positions();
        *size = cpp_positions.size();
        *positions = reinterpret_cast<chfl_vector3d*>(cpp_positions.data());
    )
}

extern "C" chfl_status chfl_frame_velocities(CHFL_FRAME* const frame, chfl_vector3d** velocities, uint64_t* size) {
    CHECK_POINTER(frame);
    CHECK_POINTER(velocities);
    CHECK_POINTER(size);
    static_assert(
        sizeof(chfl_vector3d) == sizeof(Vector3D),
        "Wrong size for chfl_vector3d. It should match Vector3D."
    );
    if (!frame->velocities()) {
        set_last_error("velocity data is not defined in this frame");
        return CHFL_MEMORY_ERROR;
    }
    CHFL_ERROR_CATCH(
        auto cpp_velocities = frame->velocities();
        *size = cpp_velocities->size();
        *velocities = reinterpret_cast<chfl_vector3d*>(cpp_velocities->data());
    )
}

extern "C" chfl_status chfl_frame_add_atom(
	CHFL_FRAME* const frame,
	const CHFL_ATOM* const atom,
	const chfl_vector3d position,
	const chfl_vector3d velocity
) {
    CHECK_POINTER(frame);
    CHECK_POINTER(atom);
    CHECK_POINTER(position);
    CHFL_ERROR_CATCH(
        if (velocity != nullptr) {
            frame->add_atom(*atom, vector3d(position), vector3d(velocity));
        } else {
            frame->add_atom(*atom, vector3d(position));
        }
    )
}

extern "C" chfl_status chfl_frame_remove(CHFL_FRAME* const frame, uint64_t i) {
    CHECK_POINTER(frame);
    CHFL_ERROR_CATCH(
        frame->remove(checked_cast(i));
    )
}

extern "C" chfl_status chfl_frame_resize(CHFL_FRAME* const frame, uint64_t size) {
    CHECK_POINTER(frame);
    CHFL_ERROR_CATCH(
        frame->resize(checked_cast(size));
    )
}

extern "C" chfl_status chfl_frame_add_velocities(CHFL_FRAME* const frame) {
    CHECK_POINTER(frame);
    CHFL_ERROR_CATCH(
        frame->add_velocities();
    )
}

extern "C" chfl_status chfl_frame_has_velocities(const CHFL_FRAME* const frame, bool* has_velocities) {
    CHECK_POINTER(frame);
    CHECK_POINTER(has_velocities);
    CHFL_ERROR_CATCH(
        *has_velocities = bool(frame->velocities());
    )
}

extern "C" chfl_status chfl_frame_set_cell(CHFL_FRAME* const frame, const CHFL_CELL* const cell) {
    CHECK_POINTER(frame);
    CHECK_POINTER(cell);
    CHFL_ERROR_CATCH(
        frame->set_cell(*cell);
    )
}

extern "C" chfl_status chfl_frame_set_topology(CHFL_FRAME* const frame, const CHFL_TOPOLOGY* const topology) {
    CHECK_POINTER(frame);
    CHECK_POINTER(topology);
    CHFL_ERROR_CATCH(
        frame->set_topology(*topology);
    )
}

extern "C" chfl_status chfl_frame_step(const CHFL_FRAME* const frame, uint64_t* step) {
    CHECK_POINTER(frame);
    CHECK_POINTER(step);
    CHFL_ERROR_CATCH(
        *step = frame->step();
    )
}

extern "C" chfl_status chfl_frame_set_step(CHFL_FRAME* const frame, uint64_t step) {
    CHECK_POINTER(frame);
    CHFL_ERROR_CATCH(
        frame->set_step(checked_cast(step));
    )
}

extern "C" chfl_status chfl_frame_guess_bonds(CHFL_FRAME* const frame) {
    CHECK_POINTER(frame);
    CHFL_ERROR_CATCH(
        frame->guess_bonds();
    )
}

extern "C" chfl_status chfl_frame_distance(const CHFL_FRAME* const frame, uint64_t i, uint64_t j, double* distance) {
    CHECK_POINTER(frame);
    CHECK_POINTER(distance);
    CHFL_ERROR_CATCH(
        *distance = frame->distance(checked_cast(i), checked_cast(j));
    )
}

extern "C" chfl_status chfl_frame_angle(const CHFL_FRAME* const frame, uint64_t i, uint64_t j, uint64_t k, double* angle) {
    CHECK_POINTER(frame);
    CHECK_POINTER(angle);
    CHFL_ERROR_CATCH(
        *angle = frame->angle(checked_cast(i), checked_cast(j), checked_cast(k));
    )
}

extern "C" chfl_status chfl_frame_dihedral(const CHFL_FRAME* const frame, uint64_t i, uint64_t j, uint64_t k, uint64_t m, double* dihedral) {
    CHECK_POINTER(frame);
    CHECK_POINTER(dihedral);
    CHFL_ERROR_CATCH(
        *dihedral = frame->dihedral(
            checked_cast(i), checked_cast(j), checked_cast(k), checked_cast(m)
        );
    )
}

extern "C" chfl_status chfl_frame_out_of_plane(const CHFL_FRAME* const frame, uint64_t i, uint64_t j, uint64_t k, uint64_t m, double* distance) {
    CHECK_POINTER(frame);
    CHECK_POINTER(distance);
    CHFL_ERROR_CATCH(
        *distance = frame->out_of_plane(
            checked_cast(i), checked_cast(j), checked_cast(k), checked_cast(m)
        );
    )
}

extern "C" chfl_status chfl_frame_properties_count(const CHFL_FRAME* const frame, uint64_t* const count) {
    CHECK_POINTER(frame);
    CHECK_POINTER(count);
    CHFL_ERROR_CATCH(
        *count = static_cast<uint64_t>(frame->properties().size());
    )
}

extern "C" chfl_status chfl_frame_list_properties(const CHFL_FRAME* const frame, const char* names[], uint64_t count) {
    CHECK_POINTER(frame);
    CHECK_POINTER(names);
    CHFL_ERROR_CATCH(
        auto& properties = frame->properties();
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

extern "C" chfl_status chfl_frame_set_property(CHFL_FRAME* const frame, const char* name, const CHFL_PROPERTY* const property) {
    CHECK_POINTER(frame);
    CHECK_POINTER(name);
    CHECK_POINTER(property);
    CHFL_ERROR_CATCH(
        frame->set(name, *property);
    )
}

extern "C" CHFL_PROPERTY* chfl_frame_get_property(const CHFL_FRAME* const frame, const char* name) {
    CHFL_PROPERTY* property = nullptr;
    CHECK_POINTER_GOTO(frame);
    CHECK_POINTER_GOTO(name);
    CHFL_ERROR_GOTO(
        auto atom_property = frame->get(name);
        if (atom_property) {
            property = shared_allocator::make_shared<Property>(*atom_property);
        } else {
            throw property_error("can not find a property named '{}' in this frame", name);
        }
    )
    return property;
error:
    chfl_free(property);
    return nullptr;
}

extern "C" chfl_status chfl_frame_add_bond(CHFL_FRAME* const frame, uint64_t i, uint64_t j) {
    CHECK_POINTER(frame);
    CHFL_ERROR_CATCH(
        frame->add_bond(checked_cast(i), checked_cast(j));
    )
}

extern "C" chfl_status chfl_frame_bond_with_order(CHFL_FRAME* const frame, uint64_t i, uint64_t j, chfl_bond_order bond_order) {
    CHECK_POINTER(frame);
    CHFL_ERROR_CATCH(
        frame->add_bond(checked_cast(i), checked_cast(j), static_cast<Bond::BondOrder>(bond_order));
    )
}

extern "C" chfl_status chfl_frame_remove_bond(CHFL_FRAME* const frame, uint64_t i, uint64_t j) {
    CHECK_POINTER(frame);
    CHFL_ERROR_CATCH(
        frame->remove_bond(checked_cast(i), checked_cast(j));
    )
}

extern "C" chfl_status chfl_frame_clear_bonds(CHFL_FRAME* const frame) {
    CHECK_POINTER(frame);
    CHFL_ERROR_CATCH(
        frame->clear_bonds();
    )
}

extern "C" chfl_status chfl_frame_add_residue(CHFL_FRAME* const frame, const CHFL_RESIDUE* const residue) {
    CHECK_POINTER(frame);
    CHECK_POINTER(residue);
    CHFL_ERROR_CATCH(
        frame->add_residue(*residue);
    )
}
