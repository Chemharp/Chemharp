// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <string>
#include <utility>

#include "chemfiles/error_fmt.hpp"
#include "chemfiles/external/optional.hpp"
#include "chemfiles/types.hpp"
#include "chemfiles/unreachable.hpp"
#include "chemfiles/warnings.hpp"

#include "chemfiles/Property.hpp"
using namespace chemfiles;

bool Property::as_bool() const {
    if (kind_ == BOOL) {
        return bool_;
    } else {
        throw property_error(
            "can not call 'as_bool' on a {} property", kind_as_string(kind_)
        );
    }
}

double Property::as_double() const {
    if (kind_ == DOUBLE) {
        return double_;
    } else {
        throw property_error(
            "can not call 'as_double' on a {} property", kind_as_string(kind_)
        );
    }
}

Vector3D Property::as_vector3d() const {
    if (kind_ == VECTOR3D) {
        return vector3d_;
    } else {
        throw property_error(
            "can not call 'as_vector3d' on a {} property", kind_as_string(kind_)
        );
    }
}

const std::string& Property::as_string() const {
    if (kind_ == STRING) {
        return string_;
    } else {
        throw property_error(
            "can not call 'as_string' on a {} property", kind_as_string(kind_)
        );
    }
}

std::string Property::kind_as_string(Kind kind) {
    switch (kind) {
    case BOOL:
        return "bool";
    case DOUBLE:
        return "double";
    case STRING:
        return "string";
    case VECTOR3D:
        return "Vector3D";
    default:
        unreachable();
    }
}


void property_map::set(std::string name, Property value) {
    // We can not move value here, because we might need it later. C++17 solves
    // this with insert_or_assign.
    auto inserted = data_.emplace(std::move(name), value);
    if (!inserted.second) {
        inserted.first->second = std::move(value);
    }
}

optional<const Property&> property_map::get(const std::string& name) const {
    auto property = data_.find(name);
    if (property != data_.end()) {
        return property->second;
    } else {
        return nullopt;
    }
}

template<Property::Kind kind>
optional<typename property_metadata<kind>::type> property_map::get(const std::string& name) const {
    auto property = this->get(name);
    if (property) {
        if (property->kind() == kind) {
            return property_metadata<kind>::extract(*property);
        } else {
            warning("",
                "expected '{}' property to be a {}, got a {} instead",
                name, Property::kind_as_string(kind), Property::kind_as_string(property->kind())
            );
            return nullopt;
        }
    } else {
        return nullopt;
    }
}

// Explicit template instantiations
template optional<bool> property_map::get<Property::BOOL>(const std::string& name) const;
template optional<double> property_map::get<Property::DOUBLE>(const std::string& name) const;
template optional<const std::string&> property_map::get<Property::STRING>(const std::string& name) const;
template optional<Vector3D> property_map::get<Property::VECTOR3D>(const std::string& name) const;
