#ifndef CHEMFILES_BINARY_FILE
#define CHEMFILES_BINARY_FILE

#include <cstddef>
#include <cstdint>

#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "chemfiles/config.h"
#include "chemfiles/File.hpp"

static_assert(sizeof(char) == sizeof(int8_t), "char must be 8-bits");

#if defined(__unix__) || (defined(__APPLE__) && defined(__MACH__))
#if (CHEMFILES_SIZEOF_VOID_P == 8)
    // use mmap in 64-bit posix
    #define CHEMFILES_BINARY_FILE_USE_MMAP 1
#else
    // otherwise fallback to an implementation using only <cstdio>
    #define CHEMFILES_BINARY_FILE_USE_MMAP 0
#endif
#else
    #define CHEMFILES_BINARY_FILE_USE_MMAP 0
#endif

namespace chemfiles {

/// A `BinaryFile` provides facilities to read/write a few primitive types
/// from/to binary (i.e. non text) files.
///
/// Depending on the file endianness, you should use one of the two subclasses
/// of this file: `BigEndianFile` or `LittleEndianFile`. All the functions
/// convert from/to the native endianess to the file endianess.
class BinaryFile: public File {
public:
    /// Open the file at the given `path` using the given `mode`
    BinaryFile(std::string path, File::Mode mode);

    /// Open the file at the given `path` using the given `mode` as a file with
    /// the current native endianess
    static std::unique_ptr<BinaryFile> open_native(std::string path, File::Mode mode);

    virtual ~BinaryFile() noexcept override;

    BinaryFile(const BinaryFile&) = delete;
    BinaryFile& operator=(const BinaryFile&) = delete;

    BinaryFile(BinaryFile&& other) noexcept : File(std::move(other)) {
        *this = std::move(other);
    }
    BinaryFile& operator=(BinaryFile&&) noexcept;

    /// Get the current position in the file
    uint64_t tell() const;

    /// Seek to the specified `position` in the file
    void seek(uint64_t position);

    /// Skip the next `count` Bytes in the file
    void skip(uint64_t count);

    /// Get the size of the file
    uint64_t file_size();

    /// Read exactly `count` char, and store them in the `data` array
    void read_char(char* data, size_t count);
    /// Read exactly as many char as fit in the pre-allocated vector
    void read_char(std::vector<char>& data) { this->read_char(data.data(), data.size()); }
    /// Read a single char value from the file
    char read_single_char() {
        char value;
        this->read_char(&value, 1);
        return value;
    }

    /// Read exactly `count` 8-bit signed integers, and store them in the `data`
    /// array
    void read_i8(int8_t* data, size_t count) {
        this->read_char(reinterpret_cast<char*>(data), count);
    }
    /// Read exactly as many 8-bit signed integers as fit in the pre-allocated vector
    void read_i8(std::vector<int8_t>& data) { this->read_i8(data.data(), data.size()); }
    /// Read a single 8-bit signed integer from the file
    int8_t read_single_i8() {
        int8_t value;
        this->read_i8(&value, 1);
        return value;
    }

    /// Read exactly `count` 8-bit unsigned integers, and store them in the `data`
    /// array
    void read_u8(uint8_t* data, size_t count) {
        this->read_char(reinterpret_cast<char*>(data), count);
    }
    /// Read exactly as many 8-bit unsigned integers as fit in the pre-allocated vector
    void read_u8(std::vector<uint8_t>& data) { this->read_u8(data.data(), data.size()); }
    /// Read a single 8-bit unsigned integer from the file
    uint8_t read_single_u8() {
        uint8_t value;
        this->read_u8(&value, 1);
        return value;
    }

    /// Read exactly `count` 16-bit signed integers, and store them in the
    /// `data` array
    virtual void read_i16(int16_t* data, size_t count) = 0;
    /// Read exactly as many 16-bit signed integers as fit in the pre-allocated vector
    void read_i16(std::vector<int16_t>& data) { this->read_i16(data.data(), data.size()); }
    /// Read a single 16-bit signed integer from the file
    int16_t read_single_i16() {
        int16_t value;
        this->read_i16(&value, 1);
        return value;
    }

    /// Read exactly `count` 16-bit unsigned integers, and store them in the
    /// `data` array
    virtual void read_u16(uint16_t* data, size_t count) = 0;
    /// Read exactly as many 16-bit unsigned integers as fit in the pre-allocated vector
    void read_u16(std::vector<uint16_t>& data) { this->read_u16(data.data(), data.size()); }
    /// Read a single 16-bit unsigned integer from the file
    uint16_t read_single_u16() {
        uint16_t value;
        this->read_u16(&value, 1);
        return value;
    }

    /// Read exactly `count` 32-bit signed integers, and store them in the
    /// `data` array
    virtual void read_i32(int32_t* data, size_t count) = 0;
    /// Read exactly as many 32-bit signed integers as fit in the pre-allocated vector
    void read_i32(std::vector<int32_t>& data) { this->read_i32(data.data(), data.size()); }
    /// Read a single 32-bit signed integer from the file
    int32_t read_single_i32() {
        int32_t value;
        this->read_i32(&value, 1);
        return value;
    }

    /// Read exactly `count` 32-bit unsigned integers, and store them in the
    /// `data` array
    virtual void read_u32(uint32_t* data, size_t count) = 0;
    /// Read exactly as many 32-bit unsigned integers as fit in the pre-allocated vector
    void read_u32(std::vector<uint32_t>& data) { this->read_u32(data.data(), data.size()); }
    /// Read a single 32-bit unsigned integer from the file
    uint32_t read_single_u32() {
        uint32_t value;
        this->read_u32(&value, 1);
        return value;
    }

    /// Read exactly `count` 64-bit signed integers, and store them in the
    /// `data` array
    virtual void read_i64(int64_t* data, size_t count) = 0;
    /// Read exactly as many 64-bit signed integers as fit in the pre-allocated vector
    void read_i64(std::vector<int64_t>& data) { this->read_i64(data.data(), data.size()); }
    /// Read a single 64-bit signed integer from the file
    int64_t read_single_i64() {
        int64_t value;
        this->read_i64(&value, 1);
        return value;
    }

    /// Read exactly `count` 64-bit unsigned integers, and store them in the
    /// `data` array
    virtual void read_u64(uint64_t* data, size_t count) = 0;
    /// Read exactly as many 64-bit unsigned integers as fit in the pre-allocated vector
    void read_u64(std::vector<uint64_t>& data) { this->read_u64(data.data(), data.size()); }
    /// Read a single 64-bit unsigned integer from the file
    uint64_t read_single_u64() {
        uint64_t value;
        this->read_u64(&value, 1);
        return value;
    }

    /// Read exactly `count` 32-bit floating point numbers, and store them in
    /// the `data` array
    virtual void read_f32(float* data, size_t count) = 0;
    /// Read exactly as many 32-bit floating point numbers as fit in the pre-allocated vector
    void read_f32(std::vector<float>& data) { this->read_f32(data.data(), data.size()); }
    /// Read a single 32-bit floating point number from the file
    float read_single_f32() {
        float value;
        this->read_f32(&value, 1);
        return value;
    }

    /// Read exactly `count` 64-bit floating point numbers, and store them in
    /// the `data` array
    virtual void read_f64(double* data, size_t count) = 0;
    /// Read exactly as many 64-bit floating point numbers as fit in the pre-allocated vector
    void read_f64(std::vector<double>& data) { this->read_f64(data.data(), data.size()); }
    /// Read a single 64-bit floating point number from the file
    double read_single_f64() {
        double value;
        this->read_f64(&value, 1);
        return value;
    }



    /// Write exactly `count` char values taken from the `data` array to the
    /// file
    void write_char(const char* data, size_t count);
    /// Write all char from the pre-allocated vector
    void write_char(const std::vector<char>& data) {
        this->write_char(data.data(), data.size());
    }
    /// Write a single char value to the file
    void write_single_char(char value) {
        this->write_char(&value, 1);
    }

    /// Write exactly `count` 8-bit signed integers taken from the `data` array
    /// to the file
    void write_i8(const int8_t* data, size_t count) {
        this->write_char(reinterpret_cast<const char*>(data), count);
    }
    /// Write all 8-bit signed integers from the pre-allocated vector
    void write_i8(const std::vector<int8_t>& data) {
        this->write_i8(data.data(), data.size());
    }
    /// Write a single 8-bit signed integer to the file
    void write_single_i8(int8_t value) {
        this->write_i8(&value, 1);
    }

    /// Write exactly `count` 8-bit unsigned integers taken from the `data`
    /// array to the file
    void write_u8(const uint8_t* data, size_t count) {
        this->write_char(reinterpret_cast<const char*>(data), count);
    }
    /// Write all 8-bit unsigned integers from the pre-allocated vector
    void write_u8(const std::vector<uint8_t>& data) {
        this->write_u8(data.data(), data.size());
    }
    /// Write a single 8-bit unsigned integer to the file
    void write_single_u8(uint8_t value) {
        this->write_u8(&value, 1);
    }

    /// Write exactly `count` 16-bit signed integers taken from the `data`
    /// array to the file
    virtual void write_i16(const int16_t* data, size_t count) = 0;
    /// Write all 16-bit signed integers from the pre-allocated vector
    void write_i16(const std::vector<int16_t>& data) {
        this->write_i16(data.data(), data.size());
    }
    /// Write a single 16-bit signed integer to the file
    void write_single_i16(int16_t value) {
        this->write_i16(&value, 1);
    }

    /// Write exactly `count` 16-bit unsigned integers taken from the `data`
    /// array to the file
    virtual void write_u16(const uint16_t* data, size_t count) = 0;
    /// Write all 16-bit unsigned integers from the pre-allocated vector
    void write_u16(const std::vector<uint16_t>& data) {
        this->write_u16(data.data(), data.size());
    }
    /// Write a single 16-bit unsigned integer to the file
    void write_single_u16(uint16_t value) {
        this->write_u16(&value, 1);
    }

    /// Write exactly `count` 32-bit signed integers taken from the `data`
    /// array to the file
    virtual void write_i32(const int32_t* data, size_t count) = 0;
    /// Write all 32-bit signed integers from the pre-allocated vector
    void write_i32(const std::vector<int32_t>& data) {
        this->write_i32(data.data(), data.size());
    }
    /// Write a single 32-bit signed integer to the file
    void write_single_i32(int32_t value) {
        this->write_i32(&value, 1);
    }

    /// Write exactly `count` 32-bit unsigned integers taken from the `data`
    /// array to the file
    virtual void write_u32(const uint32_t* data, size_t count) = 0;
    /// Write all 32-bit unsigned integers from the pre-allocated vector
    void write_u32(const std::vector<uint32_t>& data) {
        this->write_u32(data.data(), data.size());
    }
    /// Write a single 32-bit unsigned integer to the file
    void write_single_u32(uint32_t value) {
        this->write_u32(&value, 1);
    }

    /// Write exactly `count` 64-bit signed integers taken from the `data`
    /// array to the file
    virtual void write_i64(const int64_t* data, size_t count) = 0;
    /// Write all 64-bit signed integers from the pre-allocated vector
    void write_i64(const std::vector<int64_t>& data) {
        this->write_i64(data.data(), data.size());
    }
    /// Write a single 64-bit signed integer to the file
    void write_single_i64(int64_t value) {
        this->write_i64(&value, 1);
    }

    /// Write exactly `count` 64-bit unsigned integers taken from the `data`
    /// array to the file
    virtual void write_u64(const uint64_t* data, size_t count) = 0;
    /// Write all 64-bit unsigned integers from the pre-allocated vector
    void write_u64(const std::vector<uint64_t>& data) {
        this->write_u64(data.data(), data.size());
    }
    /// Write a single 64-bit unsigned integer to the file
    void write_single_u64(uint64_t value) {
        this->write_u64(&value, 1);
    }

    /// Write exactly `count` 32-bit floating point numbers taken from the
    /// `data` array to the file
    virtual void write_f32(const float* data, size_t count) = 0;
    /// Write all 32-bit floating point numbers from the pre-allocated vector
    void write_f32(const std::vector<float>& data) {
        this->write_f32(data.data(), data.size());
    }
    /// Write a single 32-bit floating point number to the file
    void write_single_f32(float value) {
        this->write_f32(&value, 1);
    }

    /// Write exactly `count` 64-bit floating point numbers taken from the
    /// `data` array to the file
    virtual void write_f64(const double* data, size_t count) = 0;
    /// Write all 64-bit floating point numbers from the pre-allocated vector
    void write_f64(const std::vector<double>& data) {
        this->write_f64(data.data(), data.size());
    }
    /// Write a single 64-bit floating point number to the file
    void write_single_f64(double value) {
        this->write_f64(&value, 1);
    }
protected:
    std::vector<char> swap_buffer_;

private:
    /// Close the file and mmap binding, and reset all member to
    /// default/moved-from values. This is used to implement both the
    /// destructor and move assignment operator
    void close_file() noexcept;

#if CHEMFILES_BINARY_FILE_USE_MMAP
    int file_descriptor_ = -1;
    char* mmap_data_ = nullptr;
    size_t file_size_ = 0;
    size_t mmap_size_ = 0;
    int mmap_prot_ = 0;
    size_t page_size_ = 0;
    uint64_t offset_ = 0;
    uint64_t total_written_size_ = 0;
#else
    FILE* file_ = nullptr;
#endif
};

class BigEndianFile: public BinaryFile {
public:
    BigEndianFile(std::string path, File::Mode mode): BinaryFile(std::move(path), mode) {}

    virtual ~BigEndianFile() noexcept override = default;

    BigEndianFile(const BigEndianFile&) = delete;
    BigEndianFile& operator=(const BigEndianFile&) = delete;

    BigEndianFile(BigEndianFile&&) noexcept = default;
    BigEndianFile& operator=(BigEndianFile&&) noexcept = default;

    void read_i16(int16_t* data, size_t count) final;
    void read_u16(uint16_t* data, size_t count) final;
    void read_i32(int32_t* data, size_t count) final;
    void read_u32(uint32_t* data, size_t count) final;
    void read_i64(int64_t* data, size_t count) final;
    void read_u64(uint64_t* data, size_t count) final;
    void read_f32(float* data, size_t count) final;
    void read_f64(double* data, size_t count) final;

    void write_i16(const int16_t* data, size_t count) final;
    void write_u16(const uint16_t* data, size_t count) final;
    void write_i32(const int32_t* data, size_t count) final;
    void write_u32(const uint32_t* data, size_t count) final;
    void write_i64(const int64_t* data, size_t count) final;
    void write_u64(const uint64_t* data, size_t count) final;
    void write_f32(const float* data, size_t count) final;
    void write_f64(const double* data, size_t count) final;
private:
    template<typename T> inline void read_as_big_endian(T* data, size_t count);
    template<typename T> inline void write_as_big_endian(const T* data, size_t count);
};

class LittleEndianFile: public BinaryFile {
public:
    LittleEndianFile(std::string path, File::Mode mode): BinaryFile(std::move(path), mode) {}

    virtual ~LittleEndianFile() noexcept override = default;

    LittleEndianFile(const LittleEndianFile&) = delete;
    LittleEndianFile& operator=(const LittleEndianFile&) = delete;

    LittleEndianFile(LittleEndianFile&&) noexcept = default;
    LittleEndianFile& operator=(LittleEndianFile&&) noexcept = default;

    void read_i16(int16_t* data, size_t count) final;
    void read_u16(uint16_t* data, size_t count) final;
    void read_i32(int32_t* data, size_t count) final;
    void read_u32(uint32_t* data, size_t count) final;
    void read_i64(int64_t* data, size_t count) final;
    void read_u64(uint64_t* data, size_t count) final;
    void read_f32(float* data, size_t count) final;
    void read_f64(double* data, size_t count) final;

    void write_i16(const int16_t* data, size_t count) final;
    void write_u16(const uint16_t* data, size_t count) final;
    void write_i32(const int32_t* data, size_t count) final;
    void write_u32(const uint32_t* data, size_t count) final;
    void write_i64(const int64_t* data, size_t count) final;
    void write_u64(const uint64_t* data, size_t count) final;
    void write_f32(const float* data, size_t count) final;
    void write_f64(const double* data, size_t count) final;
private:
    template<typename T> inline void read_as_little_endian(T* data, size_t count);
    template<typename T> inline void write_as_little_endian(const T* data, size_t count);
};

}

#endif
