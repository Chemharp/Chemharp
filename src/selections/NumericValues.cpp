// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cassert>
#include <cstring>
#include <cstdint>
#include <cstdlib>

#include <utility>

#include "chemfiles/error_fmt.hpp"
#include "chemfiles/selections/NumericValues.hpp"

using namespace chemfiles;

NumericValues::NumericValues(NumericValues&& other) noexcept: NumericValues(0) {
    *this = std::move(other);
}

NumericValues& NumericValues::operator=(NumericValues&& other) noexcept {
    std::free(heap_);

    this->value_ = other.value_;
    this->heap_ = other.heap_;

    other.heap_ = nullptr;
    return *this;
}

NumericValues::~NumericValues() {
    // NumericValues use std::malloc/std::free instead of new/delete
    // to be able to also use realloc
    std::free(heap_);
}

void NumericValues::reserve(size_t size) {
    if (heap_ == nullptr) {
        assert(size > 0);
        heap_ = static_cast<double*>(std::malloc((size + 2) * sizeof(double)));
        if (heap_ == nullptr) {
            throw error("could not allocate memory for NumericValues");
        }

        if (value_ == EMPTY_VALUE) {
            this->set_size(0);
        } else {
            heap_[2] = value_;
            this->set_size(1);
        }

        this->set_capacity(size);
    } else {
        // do not allocate if it is not needed
        if (size < this->capacity()) {
            return;
        }

        auto *heap = static_cast<double*>(std::realloc(heap_, (size + 2) * sizeof(double)));
        if (heap == nullptr) {
            throw error("could not allocate memory for NumericValues");
        }
        heap_ = heap;
        this->set_capacity(size);
    }
}

size_t NumericValues::capacity() const {
    if (heap_ == nullptr) {
        return 1;
    } else {
        uint64_t capacity;
        std::memcpy(&capacity, &heap_[1], sizeof(uint64_t));
        return static_cast<size_t>(capacity);
    }
}

/// Get the current size of this NumericValues vector
size_t NumericValues::size() const {
    if (heap_ == nullptr) {
        if (value_ == EMPTY_VALUE) {
            return 0;
        } else {
            return 1;
        }
    } else {
        uint64_t size;
        std::memcpy(&size, &heap_[0], sizeof(uint64_t));
        return static_cast<size_t>(size);
    }
}

void NumericValues::push_back(double value) {
    auto size = this->size();
    if (heap_ == nullptr) {
        if (size == 0) {
            if (value == EMPTY_VALUE) {
                throw error("invalid value +inf as first value of NumericValues");
            }
            value_ = value;
            return;
        }
        this->reserve(64);
        // reserve deals with copying from value_ to heap_
        assert(heap_[2] == value_);
    } else if (this->capacity() == size) {
        this->reserve(2 * size);
    }
    heap_[2 + size] = value;
    this->set_size(size + 1);
}

double& NumericValues::operator[](size_t i) {
    assert(i < this->size());
    if (heap_ == nullptr) {
        assert(i == 0);
        return value_;
    } else {
        return heap_[2 + i];
    }
}

double NumericValues::operator[](size_t i) const {
    return (*const_cast<NumericValues*>(this))[i];
}

const double* NumericValues::begin() const {
    if (heap_ == nullptr) {
        return &value_;
    } else {
        return heap_ + 2;
    }
}

const double* NumericValues::end() const {
    return this->begin() + this->size();
}
