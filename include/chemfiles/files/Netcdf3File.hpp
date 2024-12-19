// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_NETCDF3_FILE_HPP
#define CHEMFILES_NETCDF3_FILE_HPP

// This file contains an implementation of netcdf3 format without using the
// unidata netcdf-c library. This reduces the amount of code included, lowering
// chemfiles overall size (this is mostly important in a WASM context); and is
// significantly faster as well. Only the subset of NetCDF 3 format required for
// Amber convention is supported. This implementation is inspired by
// `scipy.io.netcdf_file`.

#include <cstdio>
#include <cstdint>

#include <functional>
#include <map>
#include <string>
#include <vector>
#include <memory>
#include <utility>

#include "chemfiles/File.hpp"
#include "chemfiles/external/optional.hpp"
#include "chemfiles/files/BinaryFile.hpp"

namespace chemfiles {
namespace netcdf3 {

namespace constants {
    // these values are taken from the main netcdf-c header, and match the
    // specification at
    // https://cdn.earthdata.nasa.gov/conduit/upload/496/ESDS-RFC-011v2.00.pdf.
    //
    // This first set of value represent the different data types the can be
    // stored in a netcdf3 file
    constexpr int32_t NC_BYTE = 1;
    constexpr int32_t NC_CHAR = 2;
    constexpr int32_t NC_SHORT = 3;
    constexpr int32_t NC_INT = 4;
    constexpr int32_t NC_FLOAT = 5;
    constexpr int32_t NC_DOUBLE = 6;

    // These values are used as marker to indicate the kind of data in the
    // following block in the file
    constexpr int32_t NC_DIMENSION = 10;
    constexpr int32_t NC_VARIABLE = 11;
    constexpr int32_t NC_ATTRIBUTE = 12;

    // these values come from the spec, and are used to fill data which was not
    // given by the user
    constexpr int8_t  NC_FILL_BYTE   = -127;
    constexpr char    NC_FILL_CHAR   = 0;
    constexpr int16_t NC_FILL_SHORT  = -32767;
    constexpr int32_t NC_FILL_INT    = -2147483647;
    constexpr float   NC_FILL_FLOAT  = 9.9692099683868690e+36f;
    constexpr double  NC_FILL_DOUBLE = 9.9692099683868690e+36;
}

class Netcdf3File;

/// Tagged enum representation for the variables or global attributes
class Value {
public:
    enum kind_t {
        BYTE,
        SHORT,
        INT,
        FLOAT,
        DOUBLE,
        STRING,
    };

    Value(char value): byte_(value), kind_(kind_t::BYTE) {}

    Value(int16_t value): short_(value), kind_(kind_t::SHORT) {}

    Value(int32_t value): int_(value), kind_(kind_t::INT) {}

    Value(float value): float_(value), kind_(kind_t::FLOAT) {}

    Value(double value): double_(value), kind_(kind_t::DOUBLE) {}

    Value(std::string value): string_(std::move(value)), kind_(kind_t::STRING) {}

    ~Value() {
        if (kind_ == kind_t::STRING) {
            using std::string;
            string_.~string();
        }
    }

    Value(Value&& other) noexcept: Value(0) {
        *this = std::move(other);
    }

    Value& operator=(Value&& other) noexcept;

    Value(const Value& other): Value(0) {
        *this = other;
    }

    Value& operator=(const Value& other);

    /// get the string stored in this `Value` or throw an error if the value
    /// does not store a string
    const std::string& as_string() const;

    /// get the 8-bit integer stored in this `Value` or throw an error if the
    /// value does not store an 8-bit integer
    int8_t as_i8() const;

    /// get the 16-bit integer stored in this `Value` or throw an error if the
    /// value does not store an 16-bit integer
    int16_t as_i16() const;

    /// get the 32-bit integer stored in this `Value` or throw an error if the
    /// value does not store an 32-bit integer
    int32_t as_i32() const;

    /// get the 32-bit floating point stored in this `Value` or throw an error
    /// if the value does not store a 32-bit floating point
    float as_f32() const;

    /// get the 64-bit floating point stored in this `Value` or throw an error
    /// if the value does not store a 64-bit floating point
    double as_f64() const;

    kind_t kind() const {
        return kind_;
    }

private:
    union {
        /// char or byte
        int8_t byte_;
        /// short
        int16_t short_;
        /// int
        int32_t int_;
        /// float
        float float_;
        /// double
        double double_;
        /// string
        std::string string_;
    };

    kind_t kind_;
};


/// A single dimension for a variable.
struct Dimension {
    Dimension(std::string name_, int32_t size_):
        name(std::move(name_)), size(size_) {}

    /// Name of the dimension
    std::string name;
    /// size of the dimension, 0 indicates the optional record (i.e. infinite)
    /// dimension
    int32_t size;

    /// Is this dimension the record (i.e. infinite) dimension?
    bool is_record() const {
        return size == 0;
    }
};

struct VariableLayout {
    /// NetCDF type identifier
    int32_t type;
    /// non record variables: the size (in bytes) of the full array, excluding padding
    /// record variables: size in bytes of a single entry, excluding padding
    int64_t size;
    /// non record variables: the size (in bytes) of the full array, including padding
    /// record variables: size in bytes of a single entry, including padding
    int64_t size_with_padding;
    /// Offset in the file of the first byte in this variable
    int64_t offset;

    /// non record variables: number of values in the full array
    /// record variables: number of values in a single entry
    size_t count() const;
    /// User-facing name for the variable type
    std::string type_name() const;
};

/// This class represent a variable in a netcdf file.
///
/// All variables have a type & shape (corresponding to a list of dimension), as
/// well as a list of attributes (key/value pairs associated with the variable).
class Variable {
public:
    /// Get all the attributes for this variable
    const std::map<std::string, Value>& attributes() const {
        return attributes_;
    }

    /// Get all the dimensions for this variable
    const std::vector<std::shared_ptr<Dimension>>& dimensions() const {
        return dimensions_;
    }

    /// Get the attribute with the given `name` for this variable, if it exists
    optional<const Value&> attribute(const std::string& name) const;

    /// Get the type of this variable as one of the values in
    /// `chemfiles::netcdf3::constants`.
    int32_t type() const {
        return layout_.type;
    }

    /// Get the shape of this variable
    std::vector<size_t> shape() const;

    /// Is this variable a record variable, i.e. one of the dimension of the
    /// variable is the record dimension?
    bool is_record() const {
        return is_record_;
    }

    /// read the content of this variable at the given `step`, filling the
    /// `data` vector. If this variable is not a record variable `step` must be
    /// 0.
    ///
    /// @throws if the vector type does not match this variable type
    template<typename T>
    void read(size_t step, std::vector<T>& data) {
        data.resize(this->layout_.count());
        this->read(step, data.data(), data.size());
    }

    /// read the content of this variable at the given `step`, writing values
    /// from `data` to `data + count`. If this variable is not a record
    /// variable `step` must be 0.
    ///
    /// @throws if `count` is too small to fit all values for this variable
    /// @throws if the pointer type does not match this variable type
    template<typename T>
    void read(size_t step, T* data, size_t count);

    /// write the content of `data` to this variable at the given `step`. If
    /// this variable is not a record variable `step` must be 0.
    ///
    /// @throws if the vector type does not match this variable type
    /// @throws if the vector size does not match this variable size
    template<typename T>
    void write(size_t step, const std::vector<T>& data) {
        this->write(step, data.data(), data.size());
    }

    /// write the values between `data` and `data + count` to this variable at
    /// the given `step`. If this variable is not a record variable `step` must
    /// be 0.
    ///
    /// @throws if the pointer type does not match this variable type
    /// @throws if count does not match this variable size
    template<typename T>
    void write(size_t step, const T* data, size_t count);

private:
    /// Write the netcdf fill/default value to this variable at the given
    /// `step`, ignoring any "fillValue" attribute. If this variable is not a
    /// record variable `step` must be 0.
    void write_fill_value(size_t step);

    friend class Netcdf3File;
    friend class Netcdf3Builder;
    Variable(
        Netcdf3File& file,
        std::vector<std::shared_ptr<Dimension>> dimensions,
        std::map<std::string, Value> attributes,
        VariableLayout layout
    );

    bool is_record_;
    std::reference_wrapper<Netcdf3File> file_;
    std::vector<std::shared_ptr<Dimension>> dimensions_;
    std::map<std::string, Value> attributes_;

    // was this variable written to at `step == file.n_records() - 1`?
    // This member is managed by both Netcdf3File (when adding a new record)
    // and the Variable itself (when writing)
    bool written_at_last_step_ = true;

    VariableLayout layout_;
};

extern template void Variable::read(size_t step, int32_t* data, size_t count);
extern template void Variable::read(size_t step, float* data, size_t count);
extern template void Variable::read(size_t step, double* data, size_t count);

extern template void Variable::write(size_t step, const char* data, size_t count);
extern template void Variable::write(size_t step, const int32_t* data, size_t count);
extern template void Variable::write(size_t step, const float* data, size_t count);
extern template void Variable::write(size_t step, const double* data, size_t count);

/// An implementation of NetCDF version 3 (or classic) binary files, using
/// 64-bit offsets for variables.
///
/// A NetCDF 3 file contains global attributes (i.e. key/values associated with
/// the file); a set of named and sized dimensions, one of which can be
/// unlimited (called the record dimension in the code); and a set of variables.
/// Each variable is linked to some dimensions, and can have its own attributes.
/// Variables which use the record dimension are stored separately, and have a
/// different API in this implementation.
///
/// NetCDF 3 files are stored using big endian, regardless of the native system
/// endianess.
class Netcdf3File: public BigEndianFile {
public:
    Netcdf3File(std::string filename, File::Mode mode);
    ~Netcdf3File() override;

    // disable moving/copying Netcdf3File since Variable instances take a
    // pointer to the file
    Netcdf3File(Netcdf3File&&) = delete;
    Netcdf3File(const Netcdf3File&) = delete;
    Netcdf3File& operator=(Netcdf3File&&) = delete;
    Netcdf3File& operator=(const Netcdf3File&) = delete;

    /// Get all the global attribute for this file
    const std::map<std::string, Value>& attributes() const {
        return attributes_;
    }

    /// Get the global attribute with the given `name` for this file, if it exists
    optional<const Value&> attribute(const std::string& name) const;

    /// Get all the variables for this file
    const std::map<std::string, Variable>& variables() const {
        return variables_;
    }

    /// Get the variable with the given name in this file if it exists
    optional<Variable&> variable(const std::string& name);

    /// Get all the dimensions defined in this file
    const std::vector<std::shared_ptr<Dimension>>& dimensions() const {
        return dimensions_;
    }

    /// Add an empty new record to this file, increasing the record dimension by
    /// one
    void add_record();

    /// get the current number of records in the file
    uint64_t n_records() const {
        return n_records_;
    }

    /// get the size in bytes of a full record entry, including all record variables
    uint64_t record_size() const {
        return record_size_;
    }

    /// Was this file already intialized (i.e. all variables/dimensions/etc.
    /// have been defined)?
    bool initialized() const {
        return initialized_;
    }

private:
    /// skip as many bytes of padding as required from the file to align the
    /// given `size` to 4-bytes
    void skip_padding(int64_t size);
    /// write as many bytes of padding as required from the file to align the
    /// given `size` to 4-bytes
    void add_padding(int64_t size);

    /// read a "Pascal" string (size + char array, no NULL terminator) from the file
    std::string read_pascal_string();
    /// write string to the file in "Pascal" format
    void write_pascal_string(const std::string& value);

    /// read attributes at the current point in file (these could be global or
    /// variable attributes)
    std::map<std::string, Value> read_attributes();
    /// read a single attribute value from the file
    Value read_attribute_value();
    /// write a single attribute value to the file
    void write_attribute_value(const Value& value);
    /// read the header for all variables
    void read_variables();

    /// current number of records in the file
    uint64_t n_records_ = 0;
    /// size in bytes of a full record entry, including all record variables
    uint64_t record_size_ = 0;

    /// list of dimensions in this file
    std::vector<std::shared_ptr<Dimension>> dimensions_;
    /// global attributes of the file
    std::map<std::string, Value> attributes_;
    /// variables in this file
    std::map<std::string, Variable> variables_;

    // was this file initialized?
    bool initialized_ = false;

    friend class Netcdf3Builder;
};

/// Definition for variables to be added to a new NetCDF 3 file
struct VariableDefinition {
    /// NetCDF type identifier for this variable
    int32_t type;
    /// list of dimensions for this variable
    std::vector<size_t> dimensions;
    /// attributes for this variable
    std::map<std::string, Value> attributes;
};

/// `Netcdf3Builder` should be used to initialize new `Netcdf3File`, setting up
/// dimensions, attributes and variables in the file.
class Netcdf3Builder {
public:
    /// Add a new dimension with the given name and size to the builder. A size
    /// of 0 indicates a record dimension
    size_t add_dimension(std::string name, size_t size);

    /// Get the list of dimensions currently defined
    const std::vector<std::shared_ptr<Dimension>>& dimensions() const {
        return dimensions_;
    }

    /// Add a new global attribute to the builder with the give `name` and `value`
    void add_attribute(std::string name, Value value);

    /// Add a new `variable` with the give `name` to the builder
    void add_variable(std::string name, VariableDefinition definition);

    /// Initialize an opened file, moving all data from this builder to the file
    void initialize(Netcdf3File* file) &&;

private:
    std::vector<std::shared_ptr<Dimension>> dimensions_;
    std::map<std::string, Value> attributes_;
    std::map<std::string, VariableDefinition> variables_;
};

} // namespace netcdf3
} // namespace chemfiles

#endif
