// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_PARSE_HPP
#define CHEMFILES_PARSE_HPP

#include <cstddef>
#include <cstdint>
#include <string>
#include <limits>
#include <string_view>
#include <type_traits>

#include "chemfiles/utils.hpp"
#include "chemfiles/error_fmt.hpp"

namespace chemfiles {

/// Convert a input to `T`, throwing a `chemfiles::Error` if the input is not
/// a valid `T`.
template<typename T> T parse(std::string_view input);

namespace detail {
    /// Helper for the static_assert below
    template<typename T>
    struct always_false { enum { value = false }; };

    template<typename Small, typename Large>
    inline Small convert_integer(Large value) {
        if (sizeof(Small) < sizeof(Large)) {
            if (value > static_cast<Large>(std::numeric_limits<Small>::max())) {
                throw error("{} is out of range for this type", value);
            }
        }
        return static_cast<Small>(value);
    }

    template<typename T>
    inline std::enable_if_t<std::is_same<T, char>::value || !std::is_integral<T>::value, T>
    parse_integer(std::string_view input) {
        (void)input;
        static_assert(
            detail::always_false<T>::value,
            "can not use parse with this type"
        );
    }

    // Conversion for all SIGNED integer type, except char
    template<typename T>
    inline std::enable_if_t<!std::is_same<T, char>::value && std::is_signed<T>::value, T>
    parse_integer(std::string_view input) {
        auto value = parse<int64_t>(input);
        return detail::convert_integer<T>(value);
    }

    // Conversion for all UNSIGNED integer type, except char
    template<typename T>
    inline std::enable_if_t<!std::is_same<T, char>::value && std::is_unsigned<T>::value, T>
    parse_integer(std::string_view input) {
        auto value = parse<uint64_t>(input);
        return detail::convert_integer<T>(value);
    }

    /// Iterator over whitespace separated values in a string
    class tokens_iterator {
    public:
        explicit tokens_iterator(std::string_view input):
            input_(input), initial_start_(input_.data()) {}

        /// Get the number of characters read from input
        size_t read_count() const {
            return static_cast<size_t>(input_.data() - initial_start_);
        }

        /// Get the next non-whitespace value. If all values have been read,
        /// this returns an empty string.
        std::string_view next() {
            auto start = input_.begin();  // NOLINT(readability-qualified-auto)
            auto end = input_.end();  // NOLINT(readability-qualified-auto)

            // skip whitespace
            while (start != end && is_ascii_whitespace(*start)) {
                start++;
            }
            input_.remove_prefix(static_cast<size_t>(start - input_.begin()));

            // Find next whitespace
            auto stop = start;  // NOLINT(readability-qualified-auto)
            while (stop != end && !is_ascii_whitespace(*stop)) {
                stop++;
            }
            auto size = static_cast<size_t>(stop - start);

            if (size == 0) {
                throw error(
                    "expected {} values, found {}",
                    count_ + 1, count_
                );
            }

            auto result = input_.substr(0, size);
            input_.remove_prefix(size);
            count_++;

            return result;
        }

    private:
        std::string_view input_;
        const char* initial_start_;
        size_t count_ = 0;
    };

    template<typename T>
    inline void scan_impl(tokens_iterator& input, T& arg) {
        arg = std::move(parse<T>(input.next()));
    }

    template<typename First, typename ...Args>
    inline void scan_impl(tokens_iterator& input, First& first, Args& ...tail) {
        first = std::move(parse<First>(input.next()));
        scan_impl(input, tail...);
    }
}

/// Read a value of type T from the `input`. This is specialized/implemented for
/// double, std::string, and all signed and unsigned integer types.
///
/// @throw chemfiles::Error if the input is empty
template<typename T>
inline T parse(std::string_view input) {
    return detail::parse_integer<T>(input);
}

/// Read a string value from the `input`. This function directly returns its
/// input.
///
/// @throw chemfiles::Error if the input is empty
template<> inline std::string parse(std::string_view input) {
    if (input.empty()) {
        throw error("tried to read a string, got an empty value");
    }
    return std::string(input);
}

/// Read double value from the `input`. This only support plain numbers (no
/// hex or octal notation), with ASCII digits (the system locale is ignored).
/// This does not support parsing NaN or infinity doubles, since they don't
/// have much interest in chemfiles;
/// Numbers should follow the `(+|-)?(\d+)?(\.\d+)?((e|E)?(+|-)?\d+)` pattern.
///
/// @throw chemfiles::Error if the input is empty, the number invalid or would
///                         overflow `double`, or if their is additional data
///                         after the value
template<> double parse(std::string_view input);

/// Read a signed 64-bit integer from the `input`. This only support plain
/// numbers (no hex or octal notation), with ASCII digits (the system locale is
/// ignored). Numbers should follow the `(+|-)?\d+` pattern.
///
/// @throw chemfiles::Error if the input is empty, the number invalid or would
///                         overflow `int64_t`, or if their is additional data
///                         after the value
template<> int64_t parse(std::string_view input);

/// Read an unsigned 64-bit integer from the `input`. This only support plain
/// numbers (no hex or octal notation), with ASCII digits (the system locale is
/// ignored). Numbers should follow the `+?\d+` pattern.
///
/// @throw chemfiles::Error if the input is empty, the number invalid or would
///                         overflow `uint64_t`, or if their is additional data
///                         after the value
template<> uint64_t parse(std::string_view input);

template<typename ...Args>
inline size_t scan(std::string_view input, Args& ...args) {
    auto iterator = detail::tokens_iterator(input);
    try {
        detail::scan_impl(iterator, args...);
    } catch (const chemfiles::Error& e) {
        throw error("error while reading '{}': {}", input, e.what());
    }
    return iterator.read_count();
}

/// Encodes an integer using the [hybrid36] encoding scheme. Returns a string
/// of `*` characters if the integer is out of range.
///
/// [hybrid36]: http://cci.lbl.gov/hybrid_36/
std::string encode_hybrid36(size_t width, int64_t value);

/// Decodes an integer using the [hybrid36] encoding scheme.
///
/// [hybrid36]: http://cci.lbl.gov/hybrid_36/
int64_t decode_hybrid36(size_t width, std::string_view input);

/// Maximum value for a width 4 number
constexpr auto MAX_HYBRID36_W4_NUMBER = 2436111;

/// Maximum value for a width 5 number
constexpr auto MAX_HYBRID36_W5_NUMBER = 87440031;

}

#endif
