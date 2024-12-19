// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_SELECTION_PARSER_HPP
#define CHEMFILES_SELECTION_PARSER_HPP

#include <cstddef>

#include <array>
#include <string>
#include <utility>
#include <vector>
#include <cassert>

#include "chemfiles/selections/lexer.hpp"
#include "chemfiles/selections/expr.hpp"

namespace chemfiles {
namespace selections {

/// Group together the actual arguments to a selector function
///
/// This is used to store `#1, #3, #2` in `angle(#1, #3, #2)`, and related
/// constructions
struct SelectionArguments {
    size_t count;
    constexpr static size_t MAX_ARGS = 4;
    std::array<SubSelection, MAX_ARGS> values;

    void add(const std::string& context, SubSelection selection);
};

/// A recursive descent parser for chemfiles selection language. This parser
/// does not handle selection context (`pairs: ...`) that should be striped
/// before using it.
class Parser {
public:
    /// Create a new parser for the given list of `tokens`
    Parser(std::vector<Token> tokens): tokens_(std::move(tokens)) {}

    ~Parser() = default;
    Parser(const Parser&) = default;
    Parser& operator=(const Parser&) = default;
    Parser(Parser&&) noexcept = default;
    Parser& operator=(Parser&&) = default;

    /// Parse the list of tokens and get the corresponding Ast.
    Ast parse();

private:
    Ast expression();
    Ast selector();
    Ast bool_selector();
    Ast string_selector();
    Ast math_selector();

    /// Parse Boolean and string properties, returning nullptr if none of these
    /// can not be parsed, so that they can be parsed as a mathematical
    /// expression later
    Ast bool_or_string_property();

    MathAst math_sum();
    MathAst math_product();
    MathAst math_power();
    MathAst math_value();
    // mathematical functions (cos, sin, ...)
    MathAst math_function(const std::string& name);
    // functions of atomic variables (distance(#1, #2), ...)
    MathAst math_var_function(const std::string& name);
    // Match multiple variables or sub-selections
    SelectionArguments arguments(const std::string& context);

    // Match an optional single variable and the surrounding parenthesis
    Variable variable();

    bool match(Token::Type type) {
        if (check(type)) {
            advance();
            return true;
        }
        return false;
    }

    bool finished() const {
        return peek().type() == Token::END;
    }

    Token peek() const {
        return tokens_[current_];
    }

    Token previous() const {
        assert(current_ > 0);
        return tokens_[current_ - 1];
    }

    bool check(Token::Type type) const {
        if (finished()) {
            return false;
        }
        return peek().type() == type;
    }

    Token advance() {
        if (!finished()) {
            current_++;
        }
        return previous();
    }

    std::vector<Token> tokens_;
    size_t current_ = 0;
};

}} // namespace chemfiles && namespace selections

#endif
