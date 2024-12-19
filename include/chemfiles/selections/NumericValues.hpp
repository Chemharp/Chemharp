// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_SELECTION_NUMERIC_VALUES_HPP
#define CHEMFILES_SELECTION_NUMERIC_VALUES_HPP

#include <cassert>
#include <cstdint>
#include <cstring>
#include <cmath>

namespace chemfiles {

/// Store multiple double values as in a std::vector<double>, while not
/// allocating heap memory for a single value.
///
/// In the heap-allocating case, all the data is behind a pointer, so that
/// sizeof(NumericValues) <= 2 * sizeof(void*), and NumericValues can be passed
/// in registers in x86_64.
///
/// The goal of this class is to try to get good performances in the common,
/// single value case, but was not benchmarked versus a simple
/// std::vector<double>.
class NumericValues {
public:
    /// Create an empty NumericValues vector
    NumericValues(): NumericValues(EMPTY_VALUE, nullptr) {}

    /// Initialize the NumericValues as stack-allocated data with the given
    /// single value
    explicit NumericValues(double value): NumericValues(value, nullptr) {}

    NumericValues(NumericValues&& other) noexcept;
    NumericValues& operator=(NumericValues&& other) noexcept;
    NumericValues(const NumericValues&) = delete;
    NumericValues& operator=(const NumericValues&) = delete;

    ~NumericValues();

    /// Get the current capacity of this NumericValues vector
    size_t capacity() const;
    /// Get the current size of this NumericValues vector
    size_t size() const;

    /// Reserve memory for `size` elements in this NumericValues vector
    void reserve(size_t size);
    /// Add the given value at the end of this NumericValues vector
    void push_back(double value);

    using iterator = const double*;

    /// Get an iterator to the beginning of this NumericValues vector
    iterator begin() const;
    /// Get an iterator to the end of this NumericValues vector
    iterator end() const;

    /// Get the element at index `i` this NumericValues vector
    double operator[](size_t i) const;
    /// Get the element at index `i` this NumericValues vector as a mutable reference
    double& operator[](size_t i);

private:
    // use `+inf` as sentinel value for len == 0
    static constexpr double EMPTY_VALUE = HUGE_VAL;

    NumericValues(double value, void* heap): value_(value), heap_(static_cast<double*>(heap)) {}

    void set_capacity(uint64_t value) {
        assert(heap_ != nullptr);
        std::memcpy(&heap_[1], &value, sizeof(uint64_t));
    }

    void set_size(uint64_t value) {
        assert(heap_ != nullptr);
        std::memcpy(&heap_[0], &value, sizeof(uint64_t));
    }

    static_assert(sizeof(double) == sizeof(uint64_t), "");

    double value_;
    // data is stored on the heap if this value is not nullptr
    // metadata about the allocation is stored at heap_[0] (size) and heap_[1]
    // (capacity) as `uint64_t`. The bytes of the uint64_t are stored inside
    // the double values.
    double* heap_;
};

} // namespace chemfiles

#endif
