// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_UTILS_HPP
#define CHEMFILES_UTILS_HPP

#include <cstddef>
#include <vector>
#include <string>
#include <algorithm>
#include <string_view>

namespace chemfiles {

/// Split `string` into components delimited by `delim`, ignoring empty
/// components.
inline std::vector<std::string_view> split(std::string_view string, char delim) {
    std::vector<std::string_view> elems;
    size_t last = 0;
    for (size_t i = 0; i<string.length(); i++) {
        if (string[i] == delim) {
            if (last != i) {
                // Don't add empty elements
                elems.push_back(string.substr(last, i - last));
            }
            last = i + 1;
        }
    }

    if (last < string.length()) {
        elems.push_back(string.substr(last, string.length() - last));
    }

    return elems;
}

inline std::vector<std::string_view> split(const char* string, char delim) {
    return split(std::string_view(string), delim);
}

// disallow temporary string
std::vector<std::string_view> split(std::string&& string, char delim) = delete;

// Check whether the given character is an ASCII whitespace
inline bool is_ascii_whitespace(char c) {
    return c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '\x0C';
}

// Check whether the given character is an ASCII lowercase letter (a-z)
inline bool is_ascii_lowercase(char c) {
    return 'a' <= c && c <= 'z';
}

// Check whether the given character is an ASCII uppercase letter (A-Z)
inline bool is_ascii_uppercase(char c) {
    return 'A' <= c && c <= 'Z';
}

// Check whether the given character is an ASCII letter (a-z, A-Z)
inline bool is_ascii_letter(char c) {
    return is_ascii_lowercase(c) || is_ascii_uppercase(c);
}

// Check whether the given character is an ASCII digit (0-9)
inline bool is_ascii_digit(char c) {
    return '0' <= c && c <= '9';
}

// Check whether the given character is an ASCII digit (0-9) or an ASCII letter
// (a-z, A-Z)
inline bool is_ascii_alphanumeric(char c) {
    return is_ascii_letter(c) || is_ascii_digit(c);
}

inline char to_ascii_lowercase(char c) {
    if (is_ascii_uppercase(c)) {
        return static_cast<char>(c - ('Z' - 'z'));
    } else {
        return c;
    }
}

inline char to_ascii_uppercase(char c) {
    if (is_ascii_lowercase(c)) {
        return static_cast<char>(c + ('Z' - 'z'));
    } else {
        return c;
    }
}

/// Remove whitespaces at the beginning and end of `string`
inline std::string_view trim(std::string_view string) {
    const auto* begin = string.data();
    const auto* end = string.data() + string.length();
    while (begin != end && is_ascii_whitespace(*begin)) {
        begin++;
    }

    if (begin != end) {
        end--;
        while (end != begin && is_ascii_whitespace(*end)) {
            end--;
        }
        end++;
    }

    return std::string_view(begin, static_cast<size_t>(end - begin));
}

inline std::string_view trim(const char* string) {
    return trim(std::string_view(string));
}

// disallow temporary string
std::string_view trim(std::string&& string) = delete;

/// Transform all characters in ASCII range in the given `string` to lower case.
///
/// Non letters and characters outside of ASCII will be left untouched
inline void to_ascii_lowercase(std::string& input) {
    std::transform(input.begin(), input.end(), input.begin(),
        [](std::string::value_type c) {
            return to_ascii_lowercase(c);
        }
    );
}

/// Transform all characters in ASCII range in the given `string` to upper case.
///
/// Non letters and characters outside of ASCII will be left untouched
inline void to_ascii_uppercase(std::string& input) {
    std::transform(input.begin(), input.end(), input.begin(),
        [](std::string::value_type c) {
            return to_ascii_uppercase(c);
        }
    );
}

/// Get the name of the computer used
std::string hostname();
/// Get the user name
std::string user_name();
/// Get the process current directory
std::string current_directory();

}

#endif
