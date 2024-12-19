// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license
#include <cassert>

#include <cstring>
#include <string>
#include <string_view>

#include <fmt/format.h>

#include "chemfiles/utils.hpp"
#include "chemfiles/error_fmt.hpp"
#include "chemfiles/external/optional.hpp"

#include "chemfiles/FormatMetadata.hpp"
using namespace chemfiles;

static void check_not_empty(const char* field_value, std::string_view field_name, const char* format_name) {
    assert(field_value != nullptr);
    if (std::string(field_value).empty()) {
        if (format_name == nullptr || field_name == "name") {
            throw format_error("the format name can not be an empty string");
        } else {
            throw format_error("the {} can not be an empty string for format '{}'", field_name, format_name);
        }
    }
}

static void check_trimmed(const char* field_value, std::string_view field_name, const char* format_name) {
    assert(field_value != nullptr);
    assert(format_name != nullptr);
    if (field_value != trim(field_value)) {
        throw format_error("the {} can not start or end with spaces for format '{}'", field_name, format_name);
    }
}

static void check_not_null(const char* field_value, std::string_view field_name, const char* format_name) {
    assert(format_name != nullptr);
    if (field_value == nullptr) {
        throw format_error("the {} can not be null for format '{}'", field_name, format_name);
    }
}


void chemfiles::FormatMetadata::validate() const {
    check_not_null(name, "name", name);
    check_not_empty(name, "name", name);
    check_trimmed(name, "name", name);

    check_not_null(description, "description", name);
    check_trimmed(description, "description", name);

    if (extension) {
        check_not_null(*extension, "extension", name);
        check_not_empty(*extension, "extension", name);
        check_trimmed(*extension, "extension", name);

        if (std::strlen(*extension) < 1 || (*extension)[0] != '.') {
            throw format_error("the extension for format '{}' must start with a dot", this->name);
        }
    }

    check_not_null(reference, "reference", name);
    check_trimmed(reference, "reference", name);

    auto reference_s = std::string(reference);
    if (!reference_s.empty()) {
        if (!(reference_s.substr(0, 7) == "http://" || reference_s.substr(0, 8) == "https://")) {
            throw format_error(
                "the reference for format '{}' must be an http link, got '{}'",
                this->name, reference_s
            );
        }
    }
}
