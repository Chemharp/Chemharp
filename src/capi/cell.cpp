// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <array>
#include <algorithm>

#include "chemfiles/capi/types.h"
#include "chemfiles/capi/misc.h"
#include "chemfiles/capi/utils.hpp"
#include "chemfiles/capi/shared_allocator.hpp"

#include "chemfiles/capi/cell.h"

#include "chemfiles/types.hpp"
#include "chemfiles/UnitCell.hpp"
#include "chemfiles/Frame.hpp"  // IWYU pragma: keep
using namespace chemfiles;

static_assert(sizeof(chfl_cellshape) == sizeof(int), "Wrong size for chfl_cellshape enum");

extern "C" CHFL_CELL* chfl_cell(const chfl_vector3d lengths, const chfl_vector3d angles) {
    CHFL_CELL* cell = nullptr;
    CHFL_ERROR_GOTO(
        if (lengths == nullptr) {
            cell = shared_allocator::make_shared<UnitCell>();
        } else if (angles == nullptr) {
            cell = shared_allocator::make_shared<UnitCell>(vector3d(lengths));
        } else {
            cell = shared_allocator::make_shared<UnitCell>(vector3d(lengths), vector3d(angles));
        }
    )
    return cell;
error:
    chfl_free(cell);
    return nullptr;
}

extern "C" CHFL_CELL* chfl_cell_from_matrix(const chfl_vector3d matrix[3]) {
    CHFL_CELL* cell = nullptr;
    CHECK_POINTER_GOTO(matrix);
    CHFL_ERROR_GOTO(
        auto cpp_matrix = Matrix3D::zero();
        std::copy(&matrix[0][0], &matrix[0][0] + 9, cpp_matrix[0].data());
        cell = shared_allocator::make_shared<UnitCell>(cpp_matrix);
    )
    return cell;
error:
    chfl_free(cell);
    return nullptr;
}


extern "C" CHFL_CELL* chfl_cell_from_frame(CHFL_FRAME* const frame) {
    CHFL_CELL* cell = nullptr;
    CHECK_POINTER_GOTO(frame);
    CHFL_ERROR_GOTO(
        cell = shared_allocator::shared_ptr<UnitCell>(frame, &frame->cell());
    )
    return cell;
error:
    chfl_free(cell);
    return nullptr;
}

extern "C" CHFL_CELL* chfl_cell_copy(const CHFL_CELL* const cell) {
    CHFL_CELL* new_cell = nullptr;
    CHFL_ERROR_GOTO(
        new_cell = shared_allocator::make_shared<UnitCell>(*cell);
    )
    return new_cell;
error:
    chfl_free(new_cell);
    return nullptr;
}

extern "C" chfl_status chfl_cell_volume(const CHFL_CELL* const cell, double* volume) {
    CHECK_POINTER(cell);
    CHECK_POINTER(volume);
    CHFL_ERROR_CATCH(
        *volume = cell->volume();
    )
}

extern "C" chfl_status chfl_cell_lengths(const CHFL_CELL* const cell, chfl_vector3d lengths) {
    CHECK_POINTER(cell);
    CHECK_POINTER(lengths);
    CHFL_ERROR_CATCH(
        auto cell_lengths = cell->lengths();
        std::copy(cell_lengths.data(), cell_lengths.data() + 3, lengths);
    )
}

extern "C" chfl_status chfl_cell_set_lengths(CHFL_CELL* const cell, const chfl_vector3d lengths) {
    CHECK_POINTER(cell);
    CHECK_POINTER(lengths);
    CHFL_ERROR_CATCH(
        cell->set_lengths(vector3d(lengths));
    )
}

extern "C" chfl_status chfl_cell_angles(const CHFL_CELL* const cell, chfl_vector3d angles) {
    CHECK_POINTER(cell);
    CHECK_POINTER(angles);
    CHFL_ERROR_CATCH(
        auto cell_angles = cell->angles();
        std::copy(cell_angles.data(), cell_angles.data() + 3, angles);
    )
}

extern "C" chfl_status chfl_cell_set_angles(CHFL_CELL* const cell, const chfl_vector3d angles) {
    CHECK_POINTER(cell);
    CHFL_ERROR_CATCH(
        cell->set_angles(vector3d(angles));
    )
}

extern "C" chfl_status chfl_cell_matrix(const CHFL_CELL* const cell, chfl_vector3d matrix[3]) {
    CHECK_POINTER(cell);
    CHECK_POINTER(matrix);
    CHFL_ERROR_CATCH(
        auto cell_matrix = cell->matrix();
        std::copy(cell_matrix[0].data(), cell_matrix[0].data() + 9, &matrix[0][0]);
    )
}

extern "C" chfl_status chfl_cell_shape(const CHFL_CELL* const cell, chfl_cellshape* const shape) {
    CHECK_POINTER(cell);
    CHECK_POINTER(shape);
    CHFL_ERROR_CATCH(
        *shape = static_cast<chfl_cellshape>(cell->shape());
    )
}

extern "C" chfl_status chfl_cell_set_shape(CHFL_CELL* const cell, chfl_cellshape shape) {
    CHECK_POINTER(cell);
    CHFL_ERROR_CATCH(
        cell->set_shape(static_cast<UnitCell::CellShape>(shape));
    )
}

extern "C" chfl_status chfl_cell_wrap(const CHFL_CELL* const cell, chfl_vector3d vector) {
    CHECK_POINTER(cell);
    CHECK_POINTER(vector);
    CHFL_ERROR_CATCH(
        auto result = cell->wrap(vector3d(vector));
        vector[0] = result[0];
        vector[1] = result[1];
        vector[2] = result[2];
    )
}
