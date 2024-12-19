// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "chemfiles/Frame.hpp"
#include "chemfiles/Selection.hpp"
#include "chemfiles/Connectivity.hpp"
#include "chemfiles/Topology.hpp"

#include "chemfiles/utils.hpp"
#include "chemfiles/error_fmt.hpp"
#include "chemfiles/unreachable.hpp"

#include "chemfiles/selections/lexer.hpp"
#include "chemfiles/selections/parser.hpp"
#include "chemfiles/selections/expr.hpp"

using namespace chemfiles;

//! Extract the context from the `string`, and put the selection string
//! (without the context) in `selection`
static Context get_context(const std::string& string, std::string& selection) {
    auto splitted = split(string, ':');
    if (splitted.size() == 1) {
        // Default to ATOM
        selection = string;
        return Context::ATOM;
    } else if (splitted.size() == 2) {
        selection = std::string(splitted[1]);
        auto context = trim(splitted[0]);
        if (context == "atoms" || context == "one") {
            return Context::ATOM;
        } else if (context == "pairs" || context == "two") {
            return Context::PAIR;
        } else if (context == "three") {
            return Context::THREE;
        } else if (context == "four") {
            return Context::FOUR;
        } else if (context == "bonds") {
            return Context::BOND;
        } else if (context == "angles") {
            return Context::ANGLE;
        } else if (context == "dihedrals") {
            return Context::DIHEDRAL;
        } else {
            throw selection_error(
                "unknown selection context '{}' in '{}'", context, string
            );
        }
    } else {
        throw selection_error(
            "can not get selection context in '{}', too many ':'", string
        );
    }
}

static unsigned max_variable(Context context) {
    switch (context) {
    case Context::ATOM:
        return 0;
    case Context::PAIR:
    case Context::BOND:
        return 1;
    case Context::ANGLE:
    case Context::THREE:
        return 2;
    case Context::DIHEDRAL:
    case Context::FOUR:
        return 3;
    default:
        unreachable();
    }
}

Selection::~Selection() = default;
Selection::Selection(Selection&&) noexcept = default;
Selection& Selection::operator=(Selection&&) noexcept = default;

Selection::Selection(std::string selection): selection_(std::move(selection)), ast_(nullptr) {
    std::string selection_string;
    context_ = get_context(selection_, selection_string);
    auto tokens = selections::Tokenizer(selection_string).tokenize();
    for (auto& token: tokens) {
        if (token.type() == selections::Token::VARIABLE) {
            if (token.variable() > max_variable(context_)) {
                throw selection_error(
                    "variable index {} is too big for the current context (should be <= {})",
                    token.variable() + 1, max_variable(context_) + 1
                );
            }
        }
    }
    ast_ = selections::Parser(tokens).parse();
    ast_->optimize();
}

size_t Selection::size() const {
    switch (context_) {
    case Context::ATOM:
        return 1;
    case Context::PAIR:
    case Context::BOND:
        return 2;
    case Context::THREE:
    case Context::ANGLE:
        return 3;
    case Context::FOUR:
    case Context::DIHEDRAL:
        return 4;
    default:
        unreachable();
    }
}

std::vector<size_t> Selection::list(const Frame& frame) const {
    if (size() != 1) {
        throw selection_error("can not call `Selection::list` on a multiple selection");
    }
    auto matches = evaluate(frame);
    auto res = std::vector<size_t>(matches.size());
    for (size_t i=0; i<matches.size(); i++) {
        res[i] = matches[i][0];
    }
    return res;
}

// Using a template to prevent putting the `is_match` function behind a pointer
template <typename match_checker>
std::vector<Match> evaluate_atoms(const Frame& frame, match_checker is_match) {
    auto matches = std::vector<Match>();
    for (size_t i=0; i<frame.size(); i++) {
        auto match = Match(i);
        if (is_match(frame, match)) {
            matches.emplace_back(match);
        }
    }
    return matches;
}

template <typename match_checker>
std::vector<Match> evaluate_pairs(const Frame& frame, match_checker is_match) {
    auto matches = std::vector<Match>();
    for (size_t i=0; i<frame.size(); i++) {
        for (size_t j=0; j<frame.size(); j++) {
            if (i == j) {continue;}
            auto match = Match(i, j);
            if (is_match(frame, match)) {
                matches.emplace_back(match);
            }
        }
    }
    return matches;
}

template <typename match_checker>
std::vector<Match> evaluate_bonds(const Frame& frame, match_checker is_match) {
    auto matches = std::vector<Match>();
    for (const auto& bond: frame.topology().bonds()) {
        auto match = Match(bond[0], bond[1]);
        if (is_match(frame, match)) {
            matches.emplace_back(match);
        } else {
            // We need to check the reverse bond (j, i); but only if
            // (i, j) is not a match. This allow to get all bonds, but
            // only once.
            match = Match(bond[1], bond[0]);
            if (is_match(frame, match)) {
                matches.emplace_back(match);
            }
        }
    }
    return matches;
}

template <typename match_checker>
std::vector<Match> evaluate_three(const Frame& frame, match_checker is_match) {
    auto matches = std::vector<Match>();
    for (size_t i=0; i<frame.size(); i++) {
        for (size_t j=0; j<frame.size(); j++) {
            if (i == j) {continue;}
            for (size_t k=0; k<frame.size(); k++) {
                if (i == k || j == k) {continue;}
                auto match = Match(i, j, k);
                if (is_match(frame, match)) {
                    matches.emplace_back(match);
                }
            }
        }
    }
    return matches;
}

template <typename match_checker>
std::vector<Match> evaluate_angles(const Frame& frame, match_checker is_match) {
    auto matches = std::vector<Match>();
    for (const auto& angle: frame.topology().angles()) {
        auto match = Match(angle[0], angle[1], angle[2]);
        if (is_match(frame, match)) {
            matches.emplace_back(match);
        } else {
            // We need to check the reverse angle (k, j, i); but only if
            // (i, j, k) is not a match. This allow to get all angles, but
            // only once.
            match = Match(angle[2], angle[1], angle[0]);
            if (is_match(frame, match)) {
                matches.emplace_back(match);
            }
        }
    }
    return matches;
}

template <typename match_checker>
std::vector<Match> evaluate_four(const Frame& frame, match_checker is_match) {
    auto matches = std::vector<Match>();
    for (size_t i=0; i<frame.size(); i++) {
        for (size_t j=0; j<frame.size(); j++) {
            if (i == j) {continue;}
            for (size_t k=0; k<frame.size(); k++) {
                if (i == k || j == k) {continue;}
                for (size_t m=0; m<frame.size(); m++) {
                    if (i == m || j == m || k == m) {continue;}
                    auto match = Match(i, j, k, m);
                    if (is_match(frame, match)) {
                        matches.emplace_back(match);
                    }
                }
            }
        }
    }
    return matches;
}

template <typename match_checker>
std::vector<Match> evaluate_dihedrals(const Frame& frame, match_checker is_match) {
    auto matches = std::vector<Match>();
    for (const auto& dihedral: frame.topology().dihedrals()) {
        auto match = Match(dihedral[0], dihedral[1], dihedral[2], dihedral[3]);
        if (is_match(frame, match)) {
            matches.emplace_back(match);
        } else {
            // We need to check the reverse dihedral (m, k, j, i); but only if
            // (i, j, k, m) is not a match. This allow to get all dihedrals, but
            // only once.
            match = Match(dihedral[3], dihedral[2], dihedral[1], dihedral[0]);
            if (is_match(frame, match)) {
                matches.emplace_back(match);
            }
        }
    }
    return matches;
}

std::vector<Match> Selection::evaluate(const Frame& frame) const {
    auto is_match = [this](const Frame& f, const Match& match) {
        return ast_->is_match(f, match);
    };

    ast_->clear();
    switch (context_) {
    case Context::ATOM:
        return evaluate_atoms(frame, is_match);
    case Context::PAIR:
        return evaluate_pairs(frame, is_match);
    case Context::BOND:
        return evaluate_bonds(frame, is_match);
    case Context::THREE:
        return evaluate_three(frame, is_match);
    case Context::ANGLE:
        return evaluate_angles(frame, is_match);
    case Context::FOUR:
        return evaluate_four(frame, is_match);
    case Context::DIHEDRAL:
        return evaluate_dihedrals(frame, is_match);
    default:
        unreachable();
    }
}
