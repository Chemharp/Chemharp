// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_PROPERTY_HPP
#define CHEMFILES_PROPERTY_HPP

#include <cstddef>
#include <string>
#include <map>
#include <utility>

#include "chemfiles/types.hpp"
#include "chemfiles/exports.h"
#include "chemfiles/unreachable.hpp"
#include "chemfiles/external/optional.hpp"

namespace chemfiles {

/// This class holds the data used in properties by the `Atom` and the `Frame`
/// classes. A property can have various types: bool, double, string or
/// `Vector3D`.
class CHFL_EXPORT Property final {
public:
    /// Possible types held by a property
    enum Kind {
        BOOL = 0,
        DOUBLE = 1,
        STRING = 2,
        VECTOR3D = 3,
    };

    /// Create a property holding a boolean value.
    /// @example{property/bool.cpp}
    Property(bool value): kind_(BOOL), bool_(value) {} // NOLINT: not all members are initialized

    /// Create a property holding a double value.
    /// @example{property/double.cpp}
    Property(double value): kind_(DOUBLE), double_(value) {} // NOLINT: not all members are initialized

    /// Create a property holding a `Vector3D` value.
    /// @example{property/vector3d.cpp}
    Property(Vector3D value): kind_(VECTOR3D), vector3d_(value) {} // NOLINT: not all members are initialized

    /// Create a property holding a string value.
    /// @example{property/string.cpp}
    Property(std::string value): kind_(STRING), string_(std::move(value)) {} // NOLINT: not all members are initialized

    // ==== The following construtors are here to help with overloading
    // ==== resolution.
    /// Create a property holding a string value from a const char*.
    /// @example{property/string.cpp}
    Property(const char* value): kind_(STRING), string_(value) {} // NOLINT: not all members are initialized

    /// Create a property holding a double value from an int.
    /// @example{property/double.cpp}
    Property(int value): kind_(DOUBLE), double_(static_cast<double>(value)) {} // NOLINT: not all members are initialized

    /// Create a property holding a double value from a long.
    /// @example{property/double.cpp}
    Property(long value): kind_(DOUBLE), double_(static_cast<double>(value)) {} // NOLINT: not all members are initialized

    /// Create a property holding a double value from a long long.
    /// @example{property/double.cpp}
    Property(long long value): kind_(DOUBLE), double_(static_cast<double>(value)) {} // NOLINT: not all members are initialized

    /// Create a property holding a double value from an unsigned.
    /// @example{property/double.cpp}
    Property(unsigned value): kind_(DOUBLE), double_(static_cast<double>(value)) {} // NOLINT: not all members are initialized

    /// Create a property holding a double value from an unsigned long.
    /// @example{property/double.cpp}
    Property(unsigned long value): kind_(DOUBLE), double_(static_cast<double>(value)) {} // NOLINT: not all members are initialized

    /// Create a property holding a double value from an unsigned long long.
    /// @example{property/double.cpp}
    Property(unsigned long long value): kind_(DOUBLE), double_(static_cast<double>(value)) {} // NOLINT: not all members are initialized

    Property(const Property& other): Property(false) {
        *this = other;
    }

    Property(Property&& other) noexcept: Property(false) {
        *this = std::move(other);
    }

    Property& operator=(const Property& other) {
        this->~Property();
        this->kind_ = other.kind_;
        switch (this->kind_) {
        case BOOL:
            new (&this->bool_) bool(other.bool_);
            break;
        case DOUBLE:
            new (&this->double_) double(other.double_);
            break;
        case STRING:
            new (&this->string_) std::string(other.string_);
            break;
        case VECTOR3D:
            new (&this->vector3d_) Vector3D(other.vector3d_);
            break;
        default:
            unreachable();
        }
        return *this;
    }

    Property& operator=(Property&& other) noexcept {
        this->~Property();
        this->kind_ = other.kind_;
        switch (this->kind_) {
        case BOOL:
            new (&this->bool_) bool(other.bool_);
            break;
        case DOUBLE:
            new (&this->double_) double(other.double_);
            break;
        case STRING:
            new (&this->string_) std::string(std::move(other.string_));
            break;
        case VECTOR3D:
            new (&this->vector3d_) Vector3D(other.vector3d_);
            break;
        default:
            unreachable();
        }
        return *this;
    }

    ~Property() {
        if (kind_ == STRING) {
            using std::string;
            this->string_.~string();
        }
    }

    /// Get the kind of property, *i.e.* the type of the held value
    ///
    /// @example{property/get_kind.cpp}
    Kind kind() const {
        return this->kind_;
    }

    /// Get the boolean value stored in this Property
    ///
    /// @example{property/as_bool.cpp}
    ///
    /// @throws PropertyError if this property does not hold a boolean value
    bool as_bool() const;

    /// Get the double value stored in this Property
    ///
    /// @example{property/as_double.cpp}
    ///
    /// @throws PropertyError if this property does not hold a double value
    double as_double() const;

    /// Get the Vector3D value stored in this Property
    ///
    /// @example{property/as_vector3d.cpp}
    ///
    /// @throws PropertyError if this property does not hold a Vector3D value
    Vector3D as_vector3d() const;

    /// Get the string value stored in this Property
    ///
    /// @example{property/as_string.cpp}
    ///
    /// @throws PropertyError if this property does not hold a string value
    const std::string& as_string() const;

    /// Get the given kind name as a string
    static std::string kind_as_string(Kind kind);

private:
    Kind kind_;
    union {
        bool bool_;
        std::string string_;
        double double_;
        Vector3D vector3d_;
    };

    friend bool operator==(const Property& lhs, const Property& rhs);
};

inline bool operator==(const Property& lhs, const Property& rhs) {
    if (lhs.kind_ != rhs.kind_) {
        return false;
    }
    switch (lhs.kind_) {
    case Property::BOOL:
        return lhs.bool_ == rhs.bool_;
    case Property::STRING:
        return lhs.string_ == rhs.string_;
    case Property::DOUBLE:
        return lhs.double_ == rhs.double_;
    case Property::VECTOR3D:
        return lhs.vector3d_ == rhs.vector3d_;
    default:
        unreachable();
    }
}

inline bool operator!=(const Property& lhs, const Property& rhs) {
    return !(lhs == rhs);
}

/// Various metadata associated with a given kind, accessible at compile time
template<Property::Kind kind>
struct property_metadata {
    /// Type of the property
    using type = void;
    /// Extract the value of the property
    static type extract(const Property&);
};

template<>
struct property_metadata<Property::BOOL> {
    using type = bool;
    static type extract(const Property& property) {
        return property.as_bool();
    }
};

template<>
struct property_metadata<Property::DOUBLE> {
    using type = double;
    static type extract(const Property& property) {
        return property.as_double();
    }
};

template<>
struct property_metadata<Property::STRING> {
    using type = const std::string&;
    static type extract(const Property& property) {
        return property.as_string();
    }
};

template<>
struct property_metadata<Property::VECTOR3D> {
    using type = Vector3D;
    static type extract(const Property& property) {
        return property.as_vector3d();
    }
};

/// A property map for inclusion in a `Frame`, an `Atom` or a `Residue`.
///
/// Properties are sorted internally, and iteration over the property will yield
/// properties in sorting order.
class CHFL_EXPORT property_map final {
public:
    using const_iterator = std::map<std::string, Property>::const_iterator;

    property_map() = default;
    property_map(property_map&&) = default;
    property_map& operator=(property_map&&) = default;
    property_map(const property_map&) = default;
    property_map& operator=(const property_map&) = default;

    /// Set an arbitrary property with the given `name` and `value`. If a
    /// property with this name already exist, it is replaced with the new
    /// value.
    void set(std::string name, Property value);

    /// Get the property with the given `name` if it exists.
    optional<const Property&> get(const std::string& name) const;

    template<Property::Kind kind>
    optional<typename property_metadata<kind>::type> get(const std::string& name) const;

    /// Get the number of properties in this property map
    size_t size() const {
        return data_.size();
    }

    /// Get an iterator to the first property in the property map
    const_iterator begin() const {
        return data_.begin();
    }

    /// Get an iterator to the end of properties
    const_iterator end() const {
        return data_.end();
    }

private:
    std::map<std::string, Property> data_;
    friend bool operator==(const property_map& lhs, const property_map& rhs);
};

inline bool operator==(const property_map& lhs, const property_map& rhs) {
    return lhs.data_ == rhs.data_;
}

inline bool operator!=(const property_map& lhs, const property_map& rhs) {
    return !(lhs == rhs);
}

// Declare instantiations of the typed `property_map::get` template
extern template CHFL_EXPORT optional<bool> property_map::get<Property::BOOL>(const std::string& name) const;
extern template CHFL_EXPORT optional<double> property_map::get<Property::DOUBLE>(const std::string& name) const;
extern template CHFL_EXPORT optional<const std::string&> property_map::get<Property::STRING>(const std::string& name) const;
extern template CHFL_EXPORT optional<Vector3D> property_map::get<Property::VECTOR3D>(const std::string& name) const;

} // namespace chemfiles

#endif
