// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_SELECTION_HPP
#define CHEMFILES_SELECTION_HPP

#include <cstddef>
#include <memory>
#include <string>
#include <vector>
#include <array>

#include "chemfiles/Error.hpp"
#include "chemfiles/exports.h"

namespace chemfiles {
class Frame;

namespace selections {
    class Selector;
    using Ast = std::unique_ptr<Selector>;
}

/// A match is a set of atomic indexes matching a given selection. The size of a
/// match depends on the associated selection, and can vary from 1 to
/// `MAX_MATCH_SIZE`.
///
/// @example{selection/match.cpp}
class CHFL_EXPORT Match final {
public:
    /// Maximal number of atoms in a match
    static constexpr size_t MAX_MATCH_SIZE = 4;

    template<typename ...Args>
    Match(Args ...args): data_({{args...}}), size_(sizeof...(args)) {
        static_assert(sizeof...(args) <= MAX_MATCH_SIZE,
        "`Match` size can not be bigger than MAX_MATCH_SIZE");
    }

    /// Get the `i`th atomic index in the match.
    ///
    /// @throw OutOfBounds if the index is bigger than the `size()` of the match
    const size_t& operator[](size_t i) const {
        if (i >= size_) {
            throw OutOfBounds("out of bounds indexing of Match");
        }
        return data_[i];
    }

    /// Size of the match, i.e. number of valid indexes in the match. This value
    /// is always equals to the size of the `Selection` this match comes from.
    size_t size() const {
        return size_;
    }

private:
    std::array<size_t, MAX_MATCH_SIZE> data_ = {{0}};
    size_t size_ = 0;
};

inline bool operator==(const Match& lhs, const Match& rhs) {
    if (lhs.size() != rhs.size()) {
        return false;
    }
    for (size_t i = 0; i<lhs.size(); i++) {
        if (lhs[i] != rhs[i]) {
            return false;
        }
    }
    return true;
}

inline bool operator!=(const Match& lhs, const Match& rhs) {
    return !(lhs == rhs);
}

/// Selection context, i.e. what we are matching
enum class Context {
    /// Matching a single atom
    ATOM,
    /// Matching a pair of atoms
    PAIR,
    /// Matching three atoms
    THREE,
    /// Matching four atoms
    FOUR,
    /// Matching a bond
    BOND,
    /// Matching an angle
    ANGLE,
    /// Matching a dihedral angle
    DIHEDRAL
};

/// This class allow to select atoms in a `Frame`, from a selection language.
///
/// The selection language is built by combining basic operations. Each basic
/// operation follows the `<selector>[(<variable>)] <operator> <value>`
/// structure, where `<operator>` is a comparison operator in `== != < <= > >=`.
///
/// @verbatim embed:rst:leading-slashes
/// Refer to the :ref:`selection-language` documentation to know the allowed
/// selectors and how to use them.
/// @endverbatim
class CHFL_EXPORT Selection final {
public:
    /// Create a selection using the given string.
    ///
    /// @throws SelectionError if there is a error in the selection string
    ///
    /// @example{selection/selection.cpp}
    explicit Selection(std::string selection);

    ~Selection();
    Selection(Selection&& other) noexcept;
    Selection& operator=(Selection&& other) noexcept;

    Selection(const Selection& other) = delete;
    Selection& operator=(const Selection& other) = delete;

    /// Evaluates the selection on a given `frame`. This function returns the
    /// list of matches in the frame for this selection.
    ///
    /// @example{selection/evaluate.cpp}
    std::vector<Match> evaluate(const Frame& frame) const;

    /// Evaluates a selection of size 1 on a given `frame`. This function
    /// returns the list of atomic indexes in the frame matching this selection.
    ///
    /// @throw SelectionError if the selection size is not 1.
    ///
    /// @example{selection/list.cpp}
    std::vector<size_t> list(const Frame& frame) const;

    /// Get the size of the selection, *i.e.* the number of atoms selected
    /// together.
    ///
    /// @example{selection/size.cpp}
    size_t size() const;

    /// Get the string used to build this selection
    ///
    /// @example{selection/string.cpp}
    std::string string() const {
        return selection_;
    }

private:
    /// Store the selection string that generated this selection
    std::string selection_;
    /// Selection context
    Context context_ = Context::ATOM;
    /// AST for evaluation of the selection
    selections::Ast ast_;
};
}

#endif
