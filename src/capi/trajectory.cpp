// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cstdint>
#include <cstring>
#include <string>

#include "chemfiles/capi/types.h"
#include "chemfiles/capi/misc.h"
#include "chemfiles/capi/utils.hpp"
#include "chemfiles/capi/shared_allocator.hpp"

#include "chemfiles/capi/trajectory.h"

#include "chemfiles/Error.hpp"
#include "chemfiles/Trajectory.hpp"
#include "chemfiles/external/span.hpp"
#include "chemfiles/external/optional.hpp"

using namespace chemfiles;

extern "C" CHFL_TRAJECTORY* chfl_trajectory_open(const char* path, char mode) {
    CHFL_TRAJECTORY* trajectory = nullptr;
    CHECK_POINTER_GOTO(path);
    CHFL_ERROR_GOTO(
        trajectory = shared_allocator::make_shared<Trajectory>(path, mode);
    )
    return trajectory;
error:
    chfl_trajectory_close(trajectory);
    return nullptr;
}

extern "C" CHFL_TRAJECTORY* chfl_trajectory_with_format(const char* path, char mode, const char* format) {
    CHFL_TRAJECTORY* trajectory = nullptr;
    CHECK_POINTER_GOTO(path);
    CHECK_POINTER_GOTO(format);
    CHFL_ERROR_GOTO(
        trajectory = shared_allocator::make_shared<Trajectory>(path, mode, format);
    )
    return trajectory;
error:
    chfl_trajectory_close(trajectory);
    return nullptr;
}

extern "C" CHFL_TRAJECTORY* chfl_trajectory_memory_reader(const char* memory, uint64_t size, const char* format) {
    CHFL_TRAJECTORY* trajectory = nullptr;
    CHECK_POINTER_GOTO(memory);
    CHECK_POINTER_GOTO(format);
    CHFL_ERROR_GOTO(
        trajectory = shared_allocator::make_shared<Trajectory>(
            Trajectory::memory_reader(memory, checked_cast(size), format)
        );
    )
    return trajectory;
error:
    chfl_trajectory_close(trajectory);
    return nullptr;
}

extern "C" CHFL_TRAJECTORY* chfl_trajectory_memory_writer(const char* format) {
    CHFL_TRAJECTORY* trajectory = nullptr;
    CHECK_POINTER_GOTO(format);
    CHFL_ERROR_GOTO(
        trajectory = shared_allocator::make_shared<Trajectory>(
            Trajectory::memory_writer(format)
        );
    )
    return trajectory;
error:
    chfl_trajectory_close(trajectory);
    return nullptr;
}

extern "C" chfl_status chfl_trajectory_path(const CHFL_TRAJECTORY* const trajectory, char* const path, uint64_t buffsize) {
    CHECK_POINTER(trajectory);
    CHECK_POINTER(path);
    CHFL_ERROR_CATCH(
        strncpy(path, trajectory->path().c_str(), checked_cast(buffsize) - 1);
        path[buffsize - 1] = '\0';
    )
}

extern "C" chfl_status chfl_trajectory_read_step(CHFL_TRAJECTORY* const trajectory, uint64_t step, CHFL_FRAME* const frame) {
    CHECK_POINTER(trajectory);
    CHECK_POINTER(frame);
    CHFL_ERROR_CATCH(
        *frame = trajectory->read_step(checked_cast(step));
    )
}

extern "C" chfl_status chfl_trajectory_read(CHFL_TRAJECTORY* const trajectory, CHFL_FRAME* const frame) {
    CHECK_POINTER(trajectory);
    CHECK_POINTER(frame);
    CHFL_ERROR_CATCH(
        *frame = trajectory->read();
    )
}

extern "C" chfl_status chfl_trajectory_write(CHFL_TRAJECTORY* const trajectory, const CHFL_FRAME* const frame) {
    CHECK_POINTER(trajectory);
    CHECK_POINTER(frame);
    CHFL_ERROR_CATCH(
         trajectory->write(*frame);
    )
}

extern "C" chfl_status chfl_trajectory_set_topology(CHFL_TRAJECTORY* const trajectory, const CHFL_TOPOLOGY* const topology) {
    CHECK_POINTER(trajectory);
    CHECK_POINTER(topology);
    CHFL_ERROR_CATCH(
        trajectory->set_topology(*topology);
    )
}

extern "C" chfl_status chfl_trajectory_topology_file(CHFL_TRAJECTORY* const trajectory, const char* path, const char* format) {
    CHECK_POINTER(trajectory);
    CHECK_POINTER(path);
    CHFL_ERROR_CATCH(
        if (format == nullptr) {
            format = "";
        }
        trajectory->set_topology(path, format);
    )
}

extern "C" chfl_status chfl_trajectory_set_cell(CHFL_TRAJECTORY* const trajectory, const CHFL_CELL* const cell) {
    CHECK_POINTER(trajectory);
    CHECK_POINTER(cell);
    CHFL_ERROR_CATCH(
        trajectory->set_cell(*cell);
    )
}

extern "C" chfl_status chfl_trajectory_nsteps(CHFL_TRAJECTORY* const trajectory, uint64_t* nsteps) {
    CHECK_POINTER(trajectory);
    CHECK_POINTER(nsteps);
    CHFL_ERROR_CATCH(
        *nsteps = trajectory->nsteps();
    )
}

extern "C" chfl_status chfl_trajectory_memory_buffer(const CHFL_TRAJECTORY* trajectory, const char** data, uint64_t* size) {
    CHECK_POINTER(trajectory);
    CHECK_POINTER(data);
    CHECK_POINTER(size);
    CHFL_ERROR_CATCH(
        auto block = trajectory->memory_buffer();
        if (!block) {
            throw Error("this trajectory was not opened to write to a memory buffer");
        }
        *data = block.value().data();
        *size = trajectory->memory_buffer().value().size();
    )
}

extern "C" void chfl_trajectory_close(const CHFL_TRAJECTORY* trajectory) {
    chfl_free(trajectory);
}
