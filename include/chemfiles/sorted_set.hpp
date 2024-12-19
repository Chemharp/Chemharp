// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_SORTED_SET_HPP
#define CHEMFILES_SORTED_SET_HPP

#include <cstddef>
#include <vector>
#include <utility>    // IWYU pragma: keep
#include <algorithm>  // IWYU pragma: keep

namespace chemfiles {

/// A set backed by a sorted vector, using binary search to insert/remove
/// values. This class follow the STL interface of `std::set`.
template<class T>
class sorted_set final: private std::vector<T> {
    using super = std::vector<T>;
public:
    using const_iterator = typename super::const_iterator;
    using iterator = typename super::const_iterator;
    using value_type = typename super::value_type;

    sorted_set(): super() {}
    sorted_set(const sorted_set&) = default;
    sorted_set& operator=(const sorted_set&) = default;
    sorted_set(sorted_set&&) = default;
    sorted_set& operator=(sorted_set&&) = default;

    using super::cbegin;
    using super::cend;
    using super::crbegin;
    using super::crend;
    const_iterator begin() const {return super::begin();}
    const_iterator end() const {return super::end();}
    const_iterator rbegin() const {return super::rbegin();}
    const_iterator rend() const {return super::rend();}

    const T& operator[](size_t i) const {return super::operator[](i);}

    using super::empty;
    using super::size;
    using super::max_size;
    using super::clear;
    using super::erase;

    std::pair<iterator, bool> insert(const value_type& value) {
        auto it = std::lower_bound(super::begin(), super::end(), value);
        if (it == super::end() || *it != value) {
            it = super::insert(it, value);
            return std::pair<iterator, bool>(it, true);
        } else {
            // Do nothing, we found the value
            return std::pair<iterator, bool>(it, false);
        }
    }

    std::pair<iterator, bool> insert(value_type&& value) {
        auto it = std::lower_bound(super::begin(), super::end(), value);
        if (it == super::end() || *it != value) {
            it = super::insert(it, std::move(value));
            return std::pair<iterator, bool>(it, true);
        } else {
            // Do nothing, we found the value
            return std::pair<iterator, bool>(it, false);
        }
    }

    template<class... Args>
    std::pair<iterator, bool> emplace(Args&&... args) {
        return insert(std::move(T(std::forward<Args>(args)...)));
    }

    iterator find(const T& value) const {
        auto it = std::lower_bound(begin(), end(), value);
        if (it == end() || *it != value) {
            return end();
        } else {
            return it;
        }
    }

    iterator erase(const_iterator pos) {
        // Converting a super::const_iterator to a super::iterator (removing the
        // const part), because GCC 4.8 do not provide the
        // std::vector<T>::erase(const_iterator) overload.
        auto it = super::begin() + (pos - super::cbegin());
        return super::erase(it);
    }

    /// Get the underlying vector data without a copy
    const super& as_vec() const {
        return *this;
    }

    /// Get the underlying vector data without a copy.
    ///
    /// WARNING: the vector is not const and can be modified, the user is
    /// supposed to ensure that the values inside stay ordered.
    super& as_mutable_vec() {
        return *this;
    }
};

} // namespace chemfiles

#endif
