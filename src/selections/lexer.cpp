// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>

#include <string>
#include <utility>
#include <vector>
#include <exception>
#include <string_view>

#include "chemfiles/parse.hpp"
#include "chemfiles/utils.hpp"
#include "chemfiles/error_fmt.hpp"
#include "chemfiles/unreachable.hpp"

#include "chemfiles/selections/lexer.hpp"

using namespace chemfiles;
using namespace chemfiles::selections;

static bool is_ident_component(char c) {
    return is_ascii_letter(c) || is_ascii_digit(c) || c == '_';
}

bool chemfiles::selections::is_ident(std::string_view string) {
    if (string.empty()) {
        return false;
    }
    if (!is_ascii_letter(string[0])) {
        return false;
    }
    for (auto c: string) {
        if (!is_ident_component(c)) {
            return false;
        }
    }
    return true;
}

std::string Token::as_str() const {
    switch (type_) {
    case Token::END:
        return "<end of selection>";
    case Token::LPAREN:
        return "(";
    case Token::RPAREN:
        return ")";
    case Token::LBRACKET:
        return "[";
    case Token::RBRACKET:
        return "]";
    case Token::COMMA:
        return ",";
    case Token::VARIABLE:
        return "#" + std::to_string(variable_ + 1);
    case Token::EQUAL:
        return "==";
    case Token::NOT_EQUAL:
        return "!=";
    case Token::LESS:
        return "<";
    case Token::LESS_EQUAL:
        return "<=";
    case Token::GREATER:
        return ">";
    case Token::GREATER_EQUAL:
        return ">=";
    case Token::PLUS:
        return "+";
    case Token::MINUS:
        return "-";
    case Token::STAR:
        return "*";
    case Token::SLASH:
        return "/";
    case Token::HAT:
        return "^";
    case Token::PERCENT:
        return "%";
    case Token::NOT:
        return "not";
    case Token::AND:
        return "and";
    case Token::OR:
        return "or";
    case Token::IDENT:
        return ident();
    case Token::STRING:
        return '"' + string() + '"';
    case Token::NUMBER:
        if (std::round(number()) == number()) {
            return std::to_string(std::lround(number()));
        } else {
            return std::to_string(number());
        }
    default:
        unreachable();
    }
}

std::vector<Token> Tokenizer::tokenize() {
    auto tokens = std::vector<Token>();
    while (!finished()) {
        if (match(is_ascii_whitespace)) {
            continue;
        } else if (match('(')) {
            tokens.emplace_back(Token(Token::LPAREN));
            continue;
        } else if (match(')')) {
            tokens.emplace_back(Token(Token::RPAREN));
            continue;
        } else if (match('[')) {
            tokens.emplace_back(Token(Token::LBRACKET));
            continue;
        } else if (match(']')) {
            tokens.emplace_back(Token(Token::RBRACKET));
            continue;
        } else if (match(',')) {
            tokens.emplace_back(Token(Token::COMMA));
            continue;
        } else if (match('=') && match('=')) {
            tokens.emplace_back(Token(Token::EQUAL));
            continue;
        } else if (match('!') && match('=')) {
            tokens.emplace_back(Token(Token::NOT_EQUAL));
            continue;
        } else if (match('<')) {
            if (match('=')) {
                tokens.emplace_back(Token(Token::LESS_EQUAL));
            } else {
                tokens.emplace_back(Token(Token::LESS));
            }
            continue;
        } else if (match('>')) {
            if (match('=')) {
                tokens.emplace_back(Token(Token::GREATER_EQUAL));
            } else {
                tokens.emplace_back(Token(Token::GREATER));
            }
            continue;
        } else if (match('+')) {
            tokens.emplace_back(Token(Token::PLUS));
            continue;
        } else if (match('-')) {
            tokens.emplace_back(Token(Token::MINUS));
            continue;
        } else if (match('*')) {
            tokens.emplace_back(Token(Token::STAR));
            continue;
        } else if (match('/')) {
            tokens.emplace_back(Token(Token::SLASH));
            continue;
        } else if (match('^')) {
            tokens.emplace_back(Token(Token::HAT));
            continue;
        } else if (match('%')) {
            tokens.emplace_back(Token(Token::PERCENT));
            continue;
        } else if (match('#')) {
            tokens.emplace_back(variable());
            continue;
        } else if (match('"')) {
            tokens.emplace_back(string());
            continue;
        } else if (check(is_ascii_letter)) {
            tokens.emplace_back(ident());
            continue;
        } else if (check(is_ascii_digit)) {
            tokens.emplace_back(number());
            continue;
        } else {
            throw selection_error("invalid character '{}' in '{}'", peek(), input_);
        }
    }
    tokens.emplace_back(Token(Token::END));
    return tokens;
}

Token Tokenizer::variable() {
    size_t start = current_;
    size_t count = 0;
    while (!finished()) {
        if (match(is_ascii_digit)) {
            count += 1;
        } else {
            break;
        }
    }

    if (count == 0) {
        throw selection_error("expected number after #");
    }

    size_t data = 0;
    auto number = input_.substr(start, count);
    try {
        data = parse<size_t>(number);
    } catch (const std::exception&) {
        throw selection_error("could not parse variable in '{}'", number);
    }

    if (data > static_cast<size_t>(UINT8_MAX)) {
        throw selection_error(
            "variable index #{} is too big (should be less than {})",
            data, UINT8_MAX
        );
    } else if (data == 0) {
        throw selection_error("invalid variable index #0");
    }
    return Token::variable(static_cast<uint8_t>(data - 1));
}

Token Tokenizer::string() {
    assert(previous() == '"');

    size_t start = current_;
    size_t count = 0;
    bool got_quote = false;

    while (!finished()) {
        if (match('"')) {
            got_quote = true;
            break;
        } else {
            advance();
            count += 1;
        }
    }

    if (!got_quote) {
        throw selection_error(
            "closing quote (\") not found in '{}'",
            input_.substr(start - 1, count + 1)
        );
    }

    return Token::string(input_.substr(start, count));
}

Token Tokenizer::ident() {
    size_t start = current_;
    size_t count = 0;
    while (!finished()) {
        if (match(is_ident_component)) {
            count += 1;
        } else {
            break;
        }
    }

    auto ident = input_.substr(start, count);
    if (ident == "or") {
        return Token(Token::OR);
    } else if (ident == "and") {
        return Token(Token::AND);
    } else if (ident == "not") {
        return Token(Token::NOT);
    }
    // Default identifier. This will be resolved during parsing phase
    return Token::ident(std::move(ident));
}

Token Tokenizer::number() {
    size_t start = current_;
    size_t count = 0;
    while (!finished()) {
        if (match(is_ascii_digit) || match('.')) {
            count += 1;
        } else {
            break;
        }
    }

    // Optional float exponent
    if (match('e') || match('E')) {
        count += 1;
        if (match('+') || match('-')) {
            count += 1;
        }
        while (!finished()) {
            if (match(is_ascii_digit)) {
                count += 1;
            } else {
                break;
            }
        }
    }

    // Require a separator between numbers and idents
    if (match(is_ident_component)) {
        count += 1;
        while (!finished()) {
            if (match(is_ident_component)) {
                count += 1;
            } else {
                break;
            }
        }
        throw selection_error("identifiers can not start with a digit: '{}'", input_.substr(start, count));
    }

    double value = 0;
    try {
        value = parse<double>(input_.substr(start, count));
    } catch (const Error& e) {
        throw SelectionError(e.what());
    }
    return Token::number(value);
}
