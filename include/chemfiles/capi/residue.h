// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_CHFL_RESIDUE_H
#define CHEMFILES_CHFL_RESIDUE_H

#include "chemfiles/capi/types.h"
#ifdef __cplusplus
extern "C" {
#endif

/// Create a new residue with the given `name` and residue identifier `resid`.
/// If the residue has no identifier, the special value of `uint64_t(-1)` should
/// be used.
///
/// The caller of this function should free the allocated memory using
/// `chfl_residue_free`.
///
/// @example{tests/capi/doc/chfl_residue/chfl_residue.c}
/// @return A pointer to the residue, or NULL in case of error.
///         You can use `chfl_last_error` to learn about the error.
CHFL_EXPORT CHFL_RESIDUE* chfl_residue(const char* name, uint64_t resid);

/// Get a copy of the residue at index `i` from a `topology`.
///
/// If `i` is bigger than the result of `chfl_topology_residues_count`, this
/// function will return `NULL`.
///
/// The residue index in the topology is not always the same as the residue
/// `id`.
///
/// The caller of this function should free the allocated memory using
/// `chfl_residue_free`.
///
/// @example{tests/capi/doc/chfl_residue/from_topology.c}
/// @return A pointer to the residue, or NULL in case of error.
///         You can use `chfl_last_error` to learn about the error.
CHFL_EXPORT CHFL_RESIDUE* chfl_residue_from_topology(
    const CHFL_TOPOLOGY* const topology, uint64_t i
);

/// Get a copy of the residue containing the atom at index `i` in the
/// `topology`.
///
/// This function will return `NULL` if the atom is not in a residue, or if the
/// index `i` is bigger than `chfl_topology_atoms_count`.
///
/// The caller of this function should free the allocated memory using
/// `chfl_residue_free`.
///
/// @example{tests/capi/doc/chfl_residue/for_atom.c}
/// @return A pointer to the residue, or NULL in case of error.
///         You can use `chfl_last_error` to learn about the error.
CHFL_EXPORT CHFL_RESIDUE* chfl_residue_for_atom(
    const CHFL_TOPOLOGY* const topology, uint64_t i
);

/// Get a copy of a `residue`.
///
/// The caller of this function should free the associated memory using
/// `chfl_residue_free`.
///
/// @example{tests/capi/doc/chfl_residue/copy.c}
/// @return A pointer to the new residue, or NULL in case of error.
///         You can use `chfl_last_error` to learn about the error.
CHFL_EXPORT CHFL_RESIDUE* chfl_residue_copy(const CHFL_RESIDUE* const residue);

/// Get the number of atoms in a `residue` in the integer pointed to by `size`.
///
/// @example{tests/capi/doc/chfl_residue/atoms_count.c}
/// @return The operation status code. You can use `chfl_last_error` to learn
///         about the error if the status code is not `CHFL_SUCCESS`.
CHFL_EXPORT chfl_status chfl_residue_atoms_count(
    const CHFL_RESIDUE* const residue, uint64_t* size
);

/// Get the identifier of a `residue` in the initial topology file in the
/// integer pointed to by `id`.
///
/// @example{tests/capi/doc/chfl_residue/id.c}
/// @return The operation status code. You can use `chfl_last_error` to learn
///         about the error if the status code is not `CHFL_SUCCESS`.
CHFL_EXPORT chfl_status chfl_residue_id(
    const CHFL_RESIDUE* const residue, uint64_t* id
);

/// Get the name of a `residue` in the string buffer `name`.
///
/// The buffer size must be passed in `buffsize`. This function will truncate
/// the residue name to fit in the buffer.
///
/// @example{tests/capi/doc/chfl_residue/name.c}
/// @return The operation status code. You can use `chfl_last_error` to learn
///         about the error if the status code is not `CHFL_SUCCESS`.
CHFL_EXPORT chfl_status chfl_residue_name(
    const CHFL_RESIDUE* const residue, char* name, uint64_t buffsize
);

/// Add the atom at index `i` in the `residue`.
///
/// @example{tests/capi/doc/chfl_residue/add_atom.c}
/// @return The operation status code. You can use `chfl_last_error` to learn
///         about the error if the status code is not `CHFL_SUCCESS`.
CHFL_EXPORT chfl_status chfl_residue_add_atom(
    CHFL_RESIDUE* const residue, uint64_t i
);

/// Check if the atom at index `i` is in the `residue`, and store the result in
/// `result`.
///
/// @example{tests/capi/doc/chfl_residue/contains.c}
/// @return The operation status code. You can use `chfl_last_error` to learn
///         about the error if the status code is not `CHFL_SUCCESS`.
CHFL_EXPORT chfl_status chfl_residue_contains(
    const CHFL_RESIDUE* const residue, uint64_t i, bool* result
);

/// Free the memory associated with a `residue`.
///
/// @example{tests/capi/doc/chfl_residue/chfl_residue.c}
/// @return `CHFL_SUCCESS`
CHFL_EXPORT chfl_status chfl_residue_free(CHFL_RESIDUE* const residue);

#ifdef __cplusplus
}
#endif

#endif
