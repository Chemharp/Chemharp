// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_CAPI_UTILS_H
#define CHEMFILES_CAPI_UTILS_H

#include <cstddef>
#include <cstdint>
#include <string>

#include <fmt/format.h>

#include "chemfiles/capi/types.h"
#include "chemfiles/warnings.hpp"  // IWYU pragma: keep
#include "chemfiles/types.hpp"

#include "chemfiles/Error.hpp"


namespace chemfiles {

void set_last_error(const std::string& message);

inline size_t checked_cast(uint64_t value) {
    // Cast SIZE_MAX to uint64_t to be able to do the comparaison, even on
    // 32-bits systems
    if (value > static_cast<uint64_t>(SIZE_MAX)) {
        throw chemfiles::Error("got a value too big to be represented by a size_t on this system");
    }
    return static_cast<size_t>(value);
}

/// Convert a C chfl_vector3d to a C++ Vector3D
inline Vector3D vector3d(const chfl_vector3d vector) {
    return {vector[0], vector[1], vector[2]};
}

#define CATCH_AND_RETURN(_exception_, _retval_)                                \
    catch (const chemfiles::_exception_& e) {                                  \
        set_last_error(e.what());                                              \
        chemfiles::send_warning(e.what());                                     \
        return _retval_;                                                       \
    }

#define CHECK_POINTER(_ptr_)                                                   \
    do {                                                                       \
        if ((_ptr_) == nullptr) {                                              \
            std::string message = fmt::format(                                 \
                "parameter '{}' cannot be NULL in {}", #_ptr_, __func__        \
            );                                                                 \
            set_last_error(message);                                           \
            chemfiles::send_warning(message.c_str());                          \
            return CHFL_MEMORY_ERROR;                                          \
        }                                                                      \
    } while (false)

#define CHECK_POINTER_GOTO(_ptr_)                                              \
    do {                                                                       \
        if ((_ptr_) == nullptr) {                                              \
            std::string message = fmt::format(                                 \
                "parameter '{}' cannot be NULL in {}", #_ptr_, __func__        \
            );                                                                 \
            set_last_error(message);                                           \
            chemfiles::send_warning(message.c_str());                          \
            goto error;  /* NOLINT: goto is OK here */                         \
        }                                                                      \
    } while (false)


/// Wrap `instructions` in a try/catch bloc automatically, and return a status
/// code
#define CHFL_ERROR_CATCH(_instructions_)                                       \
    try {                                                                      \
        _instructions_                                                         \
    }                                                                          \
    CATCH_AND_RETURN(FileError, CHFL_FILE_ERROR)                               \
    CATCH_AND_RETURN(MemoryError, CHFL_MEMORY_ERROR)                           \
    CATCH_AND_RETURN(FormatError, CHFL_FORMAT_ERROR)                           \
    CATCH_AND_RETURN(SelectionError, CHFL_SELECTION_ERROR)                     \
    CATCH_AND_RETURN(OutOfBounds, CHFL_OUT_OF_BOUNDS)                          \
    CATCH_AND_RETURN(PropertyError, CHFL_PROPERTY_ERROR)                       \
    CATCH_AND_RETURN(Error, CHFL_GENERIC_ERROR)                                \
    catch (const std::exception& e) {                                          \
        set_last_error(e.what());                                              \
        return CHFL_CXX_ERROR;                                                 \
    } catch (...) {                                                            \
        set_last_error("UNKNOWN ERROR");                                       \
        return CHFL_CXX_ERROR;                                                 \
    }                                                                          \
    return CHFL_SUCCESS;

/// Wrap `instructions` in a try/catch bloc automatically, and goto the
/// `error` label in case of error.
#define CHFL_ERROR_GOTO(_instructions_)                                        \
    try {                                                                      \
        _instructions_                                                         \
    } catch (const chemfiles::Error& e) {                                      \
        set_last_error(e.what());                                              \
        chemfiles::send_warning(e.what());                                     \
        goto error;  /* NOLINT: goto is OK here */                             \
    } catch (const std::exception& e) {                                        \
        set_last_error(e.what());                                              \
        goto error;  /* NOLINT: goto is OK here */                             \
    } catch (...) {                                                            \
        set_last_error("UNKNOWN ERROR");                                       \
        goto error;  /* NOLINT: goto is OK here */                             \
    }


} // namespace chemfiles

#endif
