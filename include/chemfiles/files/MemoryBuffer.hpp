// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_MEMORY_BUFFER_HPP
#define CHEMFILES_MEMORY_BUFFER_HPP

#include <cstddef>

#include "chemfiles/File.hpp"

namespace chemfiles {
class MemoryBuffer;

MemoryBuffer decompress_xz(const char*, size_t);
MemoryBuffer decompress_gz(const char*, size_t);
MemoryBuffer decompress_bz2(const char*, size_t);

/// A class for handling memory passed directly instead of through a file
/// handle. Unlike a `std::vector`, it does not assume ownership of the data
/// when initialized through a pointer and size.
class MemoryBuffer final {
public:
    /// Create a MemoryBuffer intended for writting with the given initial
    /// capacity
    MemoryBuffer(size_t initial);

    /// Create a MemoryBuffer intended for reading external memory
    MemoryBuffer(const char* memory, size_t size):
        // const_cast is safe here because we still encode the const-ness
        // through capacity == 0; meaning we will not try to write to this
        // pointer
        ptr_(const_cast<char*>(memory)), capacity_(0), len_(size) {}


    MemoryBuffer(const MemoryBuffer&) = delete;
    MemoryBuffer(MemoryBuffer&&) noexcept;
    MemoryBuffer& operator=(const MemoryBuffer&) = delete;
    MemoryBuffer& operator=(MemoryBuffer&&) noexcept;
    ~MemoryBuffer();

    /// Get the size of the buffer, i.e. the amount of data currently written
    /// or readable from the buffer.
    size_t size() const {
        return len_;
    }

    /// Get the capacity of the buffer, i.e. size of the current allocation
    size_t capacity() const {
        return capacity_;
    }

    const char* data() const {
        return ptr_;
    }

    void write(const char* data, size_t size);

    /// Try to decompress the content of this buffer with the given
    /// `compression` format
    void decompress(File::Compression compression);

private:
    /// Do we own the memory managed by this class?
    bool is_owned() const {
        return capacity_ != 0;
    }

    /// Reserve additional space for `extra` elements
    void reserve_extra(size_t extra);

    /// Set the length of the buffer to the given value
    void set_size(size_t len);

    /// Get a non-const pointer to the data. This is not an overload of `data`
    /// because otherwise this **private** overload would take precedence over
    /// the public one for non-const `MemoryBuffer`.
    char* data_mut() {
        return ptr_;
    }

    friend MemoryBuffer chemfiles::decompress_xz(const char*, size_t);
    friend MemoryBuffer chemfiles::decompress_gz(const char*, size_t);
    friend MemoryBuffer chemfiles::decompress_bz2(const char*, size_t);

    /// Start of the memory buffer
    char* ptr_;
    /// Size of the current allocation when writing / 0 when reading
    size_t capacity_;
    /// Amount of data written to the buffer when writing / total amount of
    /// data in the buffer when reading
    size_t len_;
};

}

#endif
