// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cassert>
#include <cstddef>
#include <cstdint>

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "chemfiles/File.hpp"
#include "chemfiles/files/Netcdf3File.hpp"

#include "chemfiles/error_fmt.hpp"
#include "chemfiles/external/optional.hpp"
#include "chemfiles/files/BinaryFile.hpp"
#include "chemfiles/unreachable.hpp"

using namespace chemfiles;
using namespace chemfiles::netcdf3;

// number of padding bytes required to align size on 4-bytes
static int64_t padding(int64_t size) {
    return (4 - (size % 4)) % 4;
}

// define some metadata for type T (name, netcdf type)
template<typename T> struct nc_type_info {};

template<> struct nc_type_info<char> {
    static constexpr int32_t nc_type = constants::NC_CHAR;
    static constexpr const char* name = "char";

    static constexpr void(Netcdf3File::*reader)(char*, size_t)
        = static_cast<void (Netcdf3File::*)(char*, size_t)>(&Netcdf3File::read_char);
    static constexpr void(Netcdf3File::*writer)(const char*, size_t)
        = static_cast<void (Netcdf3File::*)(const char*, size_t)>(&Netcdf3File::write_char);
};

template<> struct nc_type_info<int8_t> {
    static constexpr int32_t nc_type = constants::NC_BYTE;
    static constexpr const char* name = "byte/i8";

    static constexpr void(Netcdf3File::*reader)(int8_t*, size_t)
        = static_cast<void (Netcdf3File::*)(int8_t*, size_t)>(&Netcdf3File::read_i8);
    static constexpr void(Netcdf3File::*writer)(const int8_t*, size_t)
        = static_cast<void (Netcdf3File::*)(const int8_t*, size_t)>(&Netcdf3File::write_i8);
};

template<> struct nc_type_info<int16_t> {
    static constexpr int32_t nc_type = constants::NC_SHORT;
    static constexpr const char* name = "short/i16";

    static constexpr void(Netcdf3File::*reader)(int16_t*, size_t) = &Netcdf3File::read_i16;
    static constexpr void(Netcdf3File::*writer)(const int16_t*, size_t) = &Netcdf3File::write_i16;
};

template<> struct nc_type_info<int32_t> {
    static constexpr int32_t nc_type = constants::NC_INT;
    static constexpr const char* name = "int/i32";

    static constexpr void(Netcdf3File::*reader)(int32_t*, size_t) = &Netcdf3File::read_i32;
    static constexpr void(Netcdf3File::*writer)(const int32_t*, size_t) = &Netcdf3File::write_i32;
};

template<> struct nc_type_info<float> {
    static constexpr int32_t nc_type = constants::NC_FLOAT;
    static constexpr const char* name = "float/f32";

    static constexpr void(Netcdf3File::*reader)(float*, size_t) = &Netcdf3File::read_f32;
    static constexpr void(Netcdf3File::*writer)(const float*, size_t) = &Netcdf3File::write_f32;
};

template<> struct nc_type_info<double> {
    static constexpr int32_t nc_type = constants::NC_DOUBLE;
    static constexpr const char* name = "double/f64";

    static constexpr void(Netcdf3File::*reader)(double*, size_t) = &Netcdf3File::read_f64;
    static constexpr void(Netcdf3File::*writer)(const double*, size_t) = &Netcdf3File::write_f64;
};

/******************************************************************************/

Value& Value::operator=(Value&& other) noexcept {
    if (kind_ == kind_t::STRING) {
        using std::string;
        string_.~string();
    }

    this->kind_ = other.kind_;
    if (kind_ == kind_t::BYTE) {
        this->byte_ = other.byte_;
    } else if (kind_ == kind_t::SHORT) {
        this->short_ = other.short_;
    } else if (kind_ == kind_t::INT) {
        this->int_ = other.int_;
    } else if (kind_ == kind_t::FLOAT) {
        this->float_ = other.float_;
    } else if (kind_ == kind_t::DOUBLE) {
        this->double_ = other.double_;
    } else if (kind_ == kind_t::STRING) {
        // use placement new for string, assigning to a string requires a
        // valid string, but `this->string_` might not be valid yet
        new(&string_) std::string(std::move(other.string_));
    } else {
        unreachable();
    }

    return *this;
}

Value& Value::operator=(const Value& other) {
    if (kind_ == kind_t::STRING) {
        using std::string;
        string_.~string();
    }

    this->kind_ = other.kind_;
    if (kind_ == kind_t::BYTE) {
        this->byte_ = other.byte_;
    } else if (kind_ == kind_t::SHORT) {
        this->short_ = other.short_;
    } else if (kind_ == kind_t::INT) {
        this->int_ = other.int_;
    } else if (kind_ == kind_t::FLOAT) {
        this->float_ = other.float_;
    } else if (kind_ == kind_t::DOUBLE) {
        this->double_ = other.double_;
    } else if (kind_ == kind_t::STRING) {
        // use placement new for string, assigning to a string requires a
        // valid string, but `this->string_` might not be valid yet
        new(&string_) std::string(other.string_);
    } else {
        unreachable();
    }

    return *this;
}

static int32_t sizeof_nc_type(int32_t nc_type) {
    if (nc_type == constants::NC_BYTE || nc_type == constants::NC_CHAR) {
        return sizeof(char);
    } else if (nc_type == constants::NC_SHORT) {
        return sizeof(int16_t);
    } else if (nc_type == constants::NC_INT) {
        return sizeof(int32_t);
    } else if (nc_type == constants::NC_FLOAT) {
        return sizeof(float);
    } else if (nc_type == constants::NC_DOUBLE) {
        return sizeof(double);
    } else {
        unreachable();
    }
}

const std::string& Value::as_string() const {
    if (kind_ == kind_t::STRING) {
        return string_;
    } else {
        throw error("this Value does not store a string");
    }
}

int8_t Value::as_i8() const {
    if (kind_ == kind_t::BYTE) {
        return byte_;
    } else {
        throw error("this Value does not store a 8-bit integer");
    }
}

int16_t Value::as_i16() const {
    if (kind_ == kind_t::SHORT) {
        return short_;
    } else {
        throw error("this Value does not store a 16-bit integer");
    }
}

int32_t Value::as_i32() const {
    if (kind_ == kind_t::INT) {
        return int_;
    } else {
        throw error("this Value does not store a 32-bit integer");
    }
}

float Value::as_f32() const {
    if (kind_ == kind_t::FLOAT) {
        return float_;
    } else {
        throw error("this Value does not store a 32-bit floating point");
    }
}

double Value::as_f64() const {
    if (kind_ == kind_t::DOUBLE) {
        return double_;
    } else {
        throw error("this Value does not store a 64-bit floating point");
    }
}

size_t VariableLayout::count() const {
    return static_cast<size_t>(this->size) / static_cast<size_t>(sizeof_nc_type(this->type));
}

std::string VariableLayout::type_name() const {
    if (this->type == constants::NC_BYTE) {
        return nc_type_info<int8_t>::name;
    } else if (this->type == constants::NC_CHAR) {
        return nc_type_info<char>::name;
    } else if (this->type == constants::NC_SHORT) {
        return nc_type_info<int16_t>::name;
    } else if (this->type == constants::NC_INT) {
        return nc_type_info<int32_t>::name;
    } else if (this->type == constants::NC_FLOAT) {
        return nc_type_info<float>::name;
    } else if (this->type == constants::NC_DOUBLE) {
        return nc_type_info<double>::name;
    } else {
        unreachable();
    }
}

/******************************************************************************/

Variable::Variable(
    Netcdf3File& file,
    std::vector<std::shared_ptr<Dimension>> dimensions,
    std::map<std::string, Value> attributes,
    VariableLayout layout
):
    is_record_(false),
    file_(file),
    dimensions_(std::move(dimensions)),
    attributes_(std::move(attributes)),
    layout_(layout)
{
    auto size_from_file = layout_.size_with_padding;
    // the size coming from the file (in layout) is not reliable since it can
    // only go up to 2^32 - 1. Fortunately, we can compute the information
    // we need from the dimension sizes.
    layout_.size = sizeof_nc_type(layout_.type);
    for (const auto& dimension: dimensions_) {
        if (dimension->is_record()) {
            is_record_ = true;
        } else {
            layout_.size *= dimension->size;
        }
    }

    layout_.size_with_padding = layout_.size + padding(layout_.size);

    // sanity check w.r.t. the size stored in the file when possible
    if (layout_.size < INT32_MAX) {
        assert(
            layout_.size_with_padding == size_from_file
            && "declared size does not match dimension size for this variable"
        );
    }
}

optional<const Value&> Variable::attribute(const std::string& name) const {
    auto it = attributes_.find(name);
    if (it != attributes_.end()) {
        return it->second;
    } else {
        return nullopt;
    }
}

std::vector<size_t> Variable::shape() const {
    auto results = std::vector<size_t>();
    for (const auto& dimension: dimensions_) {
        if (dimension->is_record()) {
            results.emplace_back(this->file_.get().n_records());
        } else {
            results.emplace_back(dimension->size);
        }
    }
    return results;
}

template<typename T>
void Variable::read(size_t step, T* data, size_t count) {
    auto& file = file_.get();

    if (this->is_record()) {
        if (step >= file.n_records()) {
            throw file_error(
                "out of bounds: trying to read variable at step {}, "
                "but there are only {} steps in this file",
                step, file.n_records()
            );
        }
    } else if (step != 0) {
        throw file_error("can not read non-record variable at an other step than 0");
    }

    if (layout_.type != nc_type_info<T>::nc_type) {
        throw file_error(
            "internal error: the code tried to read {} data, but this "
            "variable contains {} values",
            nc_type_info<T>::name, layout_.type_name()
        );
    }

    if (count != layout_.count()) {
        throw file_error(
            "wrong array size in Variable::read: expected {}, got {}",
            layout_.count(), count
        );
    }

    auto begin = static_cast<uint64_t>(layout_.offset);
    begin += static_cast<uint64_t>(step) * file.record_size();
    file.seek(begin);
    (file.*nc_type_info<T>::reader)(data, count);
}

template void Variable::read(size_t step, int32_t* data, size_t count);
template void Variable::read(size_t step, float* data, size_t count);
template void Variable::read(size_t step, double* data, size_t count);

template<typename T>
void Variable::write(size_t step, const T* data, size_t count) {
    auto& file = file_.get();
    if (this->is_record()) {
        if (step >= file.n_records()) {
            throw file_error(
                "out of bounds: trying to write variable at step {}, "
                "but there are only {} steps in this file",
                step, file.n_records()
            );
        }
    } else if (step != 0) {
        throw file_error("can not write to non-record variable at an other step than 0");
    }

    if (layout_.type != nc_type_info<T>::nc_type) {
        throw file_error(
            "internal error: the code tried to write {} data, but this "
            "variable contains {} values",
            nc_type_info<T>::name, layout_.type_name()
        );
    }

    if (count != layout_.count()) {
        throw file_error(
            "wrong array size in Variable::read: expected {}, got {}",
            layout_.count(), count
        );
    }

    auto begin = static_cast<uint64_t>(layout_.offset);
    begin += static_cast<uint64_t>(step) * file.record_size();
    file.seek(begin);
    (file.*nc_type_info<T>::writer)(data, count);

    if (this->is_record()) {
        // fmt::print("written to {} at step {}", , step)
        if (step == file.n_records() - 1) {
            written_at_last_step_ = true;
        }
    } else {
        written_at_last_step_ = true;
    }
}

template void Variable::write(size_t step, const char* data, size_t count);
template void Variable::write(size_t step, const int32_t* data, size_t count);
template void Variable::write(size_t step, const float* data, size_t count);
template void Variable::write(size_t step, const double* data, size_t count);


void Variable::write_fill_value(size_t step) {
    if (layout_.type == constants::NC_INT) {
        auto data = std::vector<int32_t>(layout_.count(), constants::NC_FILL_INT);
        this->write(step, data);
    } else if (layout_.type == constants::NC_FLOAT) {
        auto data = std::vector<float>(layout_.count(), constants::NC_FILL_FLOAT);
        this->write(step, data);
    } else if (layout_.type == constants::NC_DOUBLE) {
        auto data = std::vector<double>(layout_.count(), constants::NC_FILL_DOUBLE);
        this->write(step, data);
    } else if (layout_.type == constants::NC_CHAR) {
        auto data = std::vector<char>(layout_.count(), constants::NC_FILL_CHAR);
        this->write(step, data);
    } else {
        throw file_error("unimplemented Variable::write_fill_value for type {}", layout_.type);
    }
}

/******************************************************************************/

Netcdf3File::Netcdf3File(std::string filename, File::Mode mode):
    BigEndianFile(std::move(filename), mode)
{
    if (mode == File::WRITE) {
        // nothing to do for now, the file will be intialized by a Netcdf3Builder
        return;
    }

    this->seek(0);
    auto buffer = std::string(3, '\0');

    try {
        this->read_char(buffer.data(), 3);
    } catch (const FileError& e) {
        if (mode == File::APPEND) {
            // An empty file was opened in append mode, it will be initialized
            // later by a Netcdf3Builder
            return;
        }
        throw e;
    }

    if (buffer != "CDF") {
        throw file_error("the file at '{}' is not a valid NetCDF 3 file", this->path());
    }

    char version = this->read_single_char();
    if (version != 2) {
        throw file_error("only 64-bit netcdf3 files are supported");
    }

    this->n_records_ = static_cast<uint64_t>(this->read_single_i32());

    // read dimensions
    auto header = this->read_single_i32();
    if (header != constants::NC_DIMENSION && header != 0) {
        throw file_error(
            "invalid netcdf3 file: expected a NC_DIMENSION header, got {}",
            header
        );
    }
    auto count = this->read_single_i32();

    for (int i=0; i<count; i++) {
        auto name = this->read_pascal_string();
        auto size = this->read_single_i32();
        dimensions_.push_back(std::make_shared<Dimension>(std::move(name), size));
    }

    attributes_ = this->read_attributes();
    this->read_variables();

    initialized_ = true;
}

Netcdf3File::~Netcdf3File() {
    if (this->mode() != File::READ) {
        // write fill values where needed
        for (auto& it: variables_) {
            auto& variable = it.second;
            if (variable.is_record()) {
                if (!variable.written_at_last_step_ && this->n_records_ > 0) {
                    variable.write_fill_value(static_cast<size_t>(this->n_records_ - 1));
                }
            } else {
                if (!variable.written_at_last_step_) {
                    variable.write_fill_value(0);
                }
            }
        }

        // update the number of records in the header, it is always at byte 4
        this->seek(4);
        this->write_single_i32(static_cast<int32_t>(this->n_records_));
    }
}

void Netcdf3File::skip_padding(int64_t size) {
    const auto count = static_cast<uint64_t>(padding(size));
    this->skip(count);
}

void Netcdf3File::add_padding(int64_t size) {
    for (int64_t i=0; i<padding(size); i++){
        this->write_single_char(0);
    }
}

std::string Netcdf3File::read_pascal_string() {
    auto size = static_cast<size_t>(this->read_single_i32());
    auto value = std::string(size, '\0');
    this->read_char(value.data(), size);
    this->skip_padding(static_cast<int64_t>(size));
    return value;
}

void Netcdf3File::write_pascal_string(const std::string& value) {
    // we don't want NULL in strings
    assert(std::string(value.c_str()) == value);

    this->write_single_i32(static_cast<int32_t>(value.size()));
    this->write_char(value.c_str(), value.size());
    this->add_padding(static_cast<int64_t>(value.size()));
}

std::map<std::string, Value> Netcdf3File::read_attributes() {
    auto header = this->read_single_i32();
    if (header != constants::NC_ATTRIBUTE && header != 0) {
        throw file_error(
            "invalid netcdf3 file: expected a NC_ATTRIBUTE header, got {}",
            header
        );
    }
    auto count = this->read_single_i32();

    auto attributes = std::map<std::string, Value>();
    for (int i=0; i<count; i++) {
        auto name = this->read_pascal_string();
        auto value = this->read_attribute_value();
        attributes.emplace(std::move(name), std::move(value));
    }

    return attributes;
}

Value Netcdf3File::read_attribute_value() {
    auto type = this->read_single_i32();
    auto count = static_cast<size_t>(this->read_single_i32());

    if (count != 1 && type != constants::NC_CHAR) {
        throw file_error("not implemented: attributes with more than one value");
    }

    size_t size = 0;
    auto value = Value(0);
    if (type == constants::NC_BYTE) {
        size = sizeof(char);
        value = Value(this->read_single_char());
    } else if (type == constants::NC_CHAR) {
        size = sizeof(char);
        auto str = std::string(static_cast<size_t>(count), '\0');
        this->read_char(str.data(), count);
        // should we remove trailing NULL if there is any?
        value = Value(std::move(str));
    } else if (type == constants::NC_SHORT) {
        size = sizeof(int16_t);
        value = Value(this->read_single_i16());
    } else if (type == constants::NC_INT) {
        size = sizeof(int32_t);
        value = Value(this->read_single_i32());
    } else if (type == constants::NC_FLOAT) {
        size = sizeof(float);
        value = Value(this->read_single_f32());
    } else if (type == constants::NC_DOUBLE) {
        size = sizeof(double);
        value = Value(this->read_single_f64());
    } else {
        throw file_error(
            "invalid netcdf3 file: attribute data type {} is unknown",
            type
        );
    }

    this->skip_padding(static_cast<int64_t>(count * size));

    return value;
}

void Netcdf3File::write_attribute_value(const Value& value) {
    size_t count = 1;
    size_t size = 0;

    auto kind = value.kind();
    if (kind == Value::BYTE) {
        size = sizeof(char);
        this->write_single_i32(constants::NC_BYTE);
        this->write_single_i32(static_cast<int32_t>(count));

        this->write_single_char(value.as_i8());
    } else if (kind == Value::SHORT) {
        size = sizeof(int16_t);
        this->write_single_i32(constants::NC_SHORT);
        this->write_single_i32(static_cast<int32_t>(count));

        this->write_single_i16(value.as_i16());
    } else if (kind == Value::INT) {
        size = sizeof(int32_t);
        this->write_single_i32(constants::NC_INT);
        this->write_single_i32(static_cast<int32_t>(count));

        this->write_single_i32(value.as_i32());
    } else if (kind == Value::FLOAT) {
        size = sizeof(float);
        this->write_single_i32(constants::NC_FLOAT);
        this->write_single_i32(static_cast<int32_t>(count));

        this->write_single_f32(value.as_f32());
    } else if (kind == Value::DOUBLE) {
        size = sizeof(double);
        this->write_single_i32(constants::NC_DOUBLE);
        this->write_single_i32(static_cast<int32_t>(count));

        this->write_single_f64(value.as_f64());
    } else if (kind == Value::STRING) {
        size = sizeof(char);
        this->write_single_i32(constants::NC_CHAR);

        const auto& str = value.as_string();
        count = str.size();

        // we don't want NULL in strings
        assert(std::string(str.c_str()) == str);
        this->write_single_i32(static_cast<int32_t>(str.size()));
        this->write_char(str.c_str(), str.size());
    } else {
        unreachable();
    }

    this->add_padding(static_cast<int64_t>(count * size));
}

void Netcdf3File::read_variables() {
    auto header = this->read_single_i32();
    if (header != constants::NC_VARIABLE && header != 0) {
        throw file_error(
            "invalid netcdf3 file: expected a NC_VARIABLE header, got {}",
            header
        );
    }
    auto count = this->read_single_i32();

    for (int var_i=0; var_i<count; var_i++) {
        auto name = this->read_pascal_string();

        auto n_dims = this->read_single_i32();
        auto dimensions = std::vector<std::shared_ptr<Dimension>>();
        for (int i=0; i<n_dims; i++) {
            auto id = static_cast<size_t>(this->read_single_i32());

            auto dimension = this->dimensions_[id];
            dimensions.emplace_back(dimension);
        }

        auto attributes = this->read_attributes();

        auto type = this->read_single_i32();
        auto size_with_padding = this->read_single_i32();
        // this is where the 64-bit offset changes something to the format
        // auto offset = this->read_single_i32();
        auto offset = this->read_single_i64();

        variables_.emplace(std::move(name), Variable(
            *this,
            std::move(dimensions),
            std::move(attributes),
            VariableLayout { type, /*size*/ 0, size_with_padding, offset }
        ));
    }

    // now that we know all the variables, compute the size of a single record
    // entry, which includes all record variables
    this->record_size_ = 0;
    for (const auto& it: variables_) {
        const auto& variable = it.second;

        if (variable.is_record()) {
            this->record_size_ += static_cast<uint64_t>(variable.layout_.size);
            this->record_size_ += static_cast<uint64_t>(padding(variable.layout_.size));
        }
    }
}

optional<const Value&> Netcdf3File::attribute(const std::string& name) const {
    auto it = attributes_.find(name);
    if (it != attributes_.end()) {
        return it->second;
    } else {
        return nullopt;
    }
}

optional<Variable&> Netcdf3File::variable(const std::string& name) {
    auto it = variables_.find(name);
    if (it != variables_.end()) {
        return it->second;
    } else {
        return nullopt;
    }
}

void Netcdf3File::add_record() {
    if (this->mode() != File::WRITE && this->mode() != File::APPEND) {
        throw file_error("can not add a record to a file opened in read-only mode");
    }

    this->n_records_ += 1;

    // write fill value in record variables for last step if needed & mark them
    // as not written to for this step
    for (auto& it: variables_) {
        auto& variable = it.second;
        if (variable.is_record()) {
            if (!variable.written_at_last_step_ && this->n_records_ > 1) {
                variable.write_fill_value(static_cast<size_t>(this->n_records_ - 2));
            }
            variable.written_at_last_step_ = false;
        }
    }
}

/******************************************************************************/

size_t Netcdf3Builder::add_dimension(std::string name, size_t size) {
    auto dimension = std::make_shared<Dimension>(std::move(name), size);
    if (dimension->size < 0) {
        throw file_error(
            "dimension size must be positive, got {} for '{}'",
            dimension->size, dimension->name
        );
    }

    for (const auto& existing: dimensions_) {
        if (existing->is_record() && dimension->is_record()) {
            throw file_error(
                "only one dimension can be a record dimension, already got {}",
                existing->name
            );
        }

        if (existing->name == dimension->name) {
            throw file_error(
                "can not add a dimension named '{}', already got one",
                existing->name
            );
        }

        // there is no validation that the dimension name is valid according
        // to netcdf3 spec
    }

    dimensions_.emplace_back(std::move(dimension));
    return dimensions_.size() - 1;
}

void Netcdf3Builder::add_attribute(std::string name, Value value) {
    if (attributes_.find(name) != attributes_.end()) {
        throw file_error("attribute '{}' already exists in this file", name);
    }
    attributes_.emplace(std::move(name), std::move(value));
}

void Netcdf3Builder::add_variable(std::string name, VariableDefinition definition) {
    if (variables_.find(name) != variables_.end()) {
        throw file_error("variable '{}' already exists in this file", name);
    }

    for (auto id: definition.dimensions) {
        if (id >= this->dimensions_.size()) {
            throw file_error("invalid dimension id for variable '{}'", name);
        }
    }

    if (definition.type != constants::NC_CHAR
        && definition.type != constants::NC_BYTE
        && definition.type != constants::NC_SHORT
        && definition.type != constants::NC_INT
        && definition.type != constants::NC_FLOAT
        && definition.type != constants::NC_DOUBLE
    ) {
        throw file_error("invalid type for variable '{}'", name);
    }

    variables_.emplace(std::move(name), std::move(definition));
}

void Netcdf3Builder::initialize(Netcdf3File* file) && {
    if (file->initialized_) {
        throw file_error("can not initialize a file twice");
    }
    file->initialized_ = true;

    file->seek(0);

    file->write_char("CDF", 3);
    file->write_single_char(2); // version == 2 for 64-bit offset
    file->write_single_i32(0);  // record length

    // write dimensions
    file->write_single_i32(constants::NC_DIMENSION);
    file->write_single_i32(static_cast<int32_t>(dimensions_.size()));
    for (const auto& dimension: dimensions_) {
        file->write_pascal_string(dimension->name);
        file->write_single_i32(dimension->size);
    }

    // write attributes
    file->write_single_i32(constants::NC_ATTRIBUTE);
    file->write_single_i32(static_cast<int32_t>(attributes_.size()));
    for (const auto& it: attributes_) {
        file->write_pascal_string(it.first);
        file->write_attribute_value(it.second);
    }

    auto variables = std::map<std::string, Variable>();
    // positions where the variables offset should be written
    auto offset_positions = std::map<std::string, uint64_t>();

    // write variables
    file->write_single_i32(constants::NC_VARIABLE);
    file->write_single_i32(static_cast<int32_t>(variables_.size()));
    for (auto it: std::move(variables_)) {
        auto name = it.first;
        auto variable = std::move(it.second);

        file->write_pascal_string(name);
        int64_t size_with_padding = sizeof_nc_type(variable.type);

        // write variable dimensions
        auto dimensions = std::vector<std::shared_ptr<Dimension>>();
        file->write_single_i32(static_cast<int32_t>(variable.dimensions.size()));
        for (const auto& id: variable.dimensions) {
            file->write_single_i32(static_cast<int32_t>(id));

            auto dimension = dimensions_[id];
            if (!dimension->is_record()) {
                size_with_padding *= dimension->size;
            }

            dimensions.emplace_back(dimension);
        }

        // add padding to size
        size_with_padding += padding(size_with_padding);

        // write variable attributes
        file->write_single_i32(constants::NC_ATTRIBUTE);
        file->write_single_i32(static_cast<int32_t>(variable.attributes.size()));
        for (const auto& attribute: variable.attributes) {
            file->write_pascal_string(attribute.first);
            file->write_attribute_value(attribute.second);
        }

        file->write_single_i32(variable.type);
        file->write_single_i32(static_cast<int32_t>(size_with_padding));

        // the actual offset will be set below, once all variable metadata have
        // been added to the file
        int64_t offset = -1;
        offset_positions.emplace(name, file->tell());
        file->write_single_i64(0);

        variables.emplace(std::move(name), Variable(
            *file,
            std::move(dimensions),
            std::move(variable.attributes),
            VariableLayout { variable.type, /*size*/ 0, size_with_padding, offset }
        ));
    }

    // now that we know the size taken by the file metadata, we can compute the
    // offset for all variables, starting by the non-record variables
    uint64_t offset = file->tell();
    for (auto& it: variables) {
        if (it.second.is_record()) {
            continue;
        }

        auto& layout = it.second.layout_;
        // write offset to the variable
        layout.offset = static_cast<int64_t>(offset);
        // and to the file
        file->seek(offset_positions.at(it.first));
        file->write_single_i64(static_cast<int64_t>(offset));

        auto delta = static_cast<size_t>(layout.size) + static_cast<size_t>(padding(layout.size));
        offset += delta;
    }

    uint64_t record_size = 0;
    for (auto& it: variables) {
        if (!it.second.is_record()) {
            continue;
        }

        auto& layout = it.second.layout_;
        // write offset to the variable
        layout.offset = static_cast<int64_t>(offset);
        // and to the file
        file->seek(offset_positions.at(it.first));
        file->write_single_i64(static_cast<int64_t>(offset));

        auto delta = static_cast<uint64_t>(layout.size) + static_cast<uint64_t>(padding(layout.size));
        offset += delta;
        record_size += delta;
    }

    // fill up the non-record variable with fill values
    for (auto& it: variables) {
        if (!it.second.is_record()) {
            it.second.write_fill_value(0);
        }
    }

    // move all the data to the file
    file->record_size_ = record_size;
    file->n_records_ = 0;
    file->attributes_ = std::move(this->attributes_);
    file->dimensions_ = std::move(this->dimensions_);
    file->variables_ = std::move(variables);
}
