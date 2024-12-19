// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cstdlib>
#include <cstring>
#include <utility>
#include <algorithm>

#include "chemfiles/error_fmt.hpp"
#include "chemfiles/unreachable.hpp"

#include "chemfiles/File.hpp"
#include "chemfiles/files/GzFile.hpp"
#include "chemfiles/files/XzFile.hpp"
#include "chemfiles/files/Bz2File.hpp"

#include "chemfiles/files/MemoryBuffer.hpp"

using namespace chemfiles;

MemoryBuffer::MemoryBuffer(size_t initial):
    ptr_(nullptr),
    capacity_(0),
    len_(0)
{
    if (initial == 0) {
        throw file_error("invalid initial size of 0 for MemoryBuffer");
    }
    /// use malloc/free instead of new/delete to also have access to calloc &
    /// realloc
    ptr_ = static_cast<char*>(std::calloc(initial, sizeof(char)));
    if (ptr_ == nullptr) {
        throw file_error("failed to allocate memory for MemoryBuffer");
    }
    capacity_ = initial;
}


MemoryBuffer::MemoryBuffer(MemoryBuffer&& other) noexcept: MemoryBuffer(nullptr, 0) {
    *this = std::move(other);
}

MemoryBuffer& MemoryBuffer::operator=(MemoryBuffer&& other) noexcept {
    if (this->is_owned()) {
        std::free(this->ptr_);
    }

    this->ptr_ = other.ptr_;
    this->len_ = other.len_;
    this->capacity_ = other.capacity_;

    other.ptr_ = nullptr;
    other.len_ = 0;
    other.capacity_ = 0;
    return *this;
}

MemoryBuffer::~MemoryBuffer() {
    if (this->is_owned()) {
        std::free(ptr_);
    }
}

void MemoryBuffer::write(const char* data, size_t size) {
    if (!this->is_owned()) {
        throw file_error("can not write to read-only MemoryBuffer");
    }

    // +1 to always ensure there is a NULL byte at the end of the buffer
    if (size > (capacity_ - (len_ + 1))) {
        auto extra = capacity_;
        while (size > (capacity_ + extra - (len_ + 1))) {
            extra *= 2;
        }
        this->reserve_extra(extra);
    }

    auto* start = ptr_ + len_;
    std::copy(data, data + size, start);
    len_ += size;
}

void MemoryBuffer::reserve_extra(size_t extra) {
    if (!this->is_owned()) {
        throw file_error("can not reserve additional memory for a read-only MemoryBuffer");
    }

    ptr_ = static_cast<char*>(std::realloc(ptr_, capacity_ + extra));
    if (ptr_ == nullptr) {
        std::free(ptr_);
        throw file_error("failed to allocate memory for MemoryBuffer");
    }
    // set newly allocated memory to 0
    std::memset(ptr_ + capacity_, 0, extra);
    capacity_ += extra;
}

void MemoryBuffer::set_size(size_t len) {
    if (!this->is_owned()) {
        throw file_error("can not set the size of a read-only MemoryBuffer");
    }

    if (len >= capacity_) {
        throw file_error("can not set the size of a MemoryBuffer larger than the capacity");
    }

    len_ = len;
}

void MemoryBuffer::decompress(File::Compression compression) {
    switch(compression) {
    case File::GZIP:
        *this = decompress_gz(this->data(), this->size());
        break;
    case File::LZMA:
        *this = decompress_xz(this->data(), this->size());
        break;
    case File::BZIP2:
        *this = decompress_bz2(this->data(), this->size());
        break;
    case File::DEFAULT:
        // nothing to do
        break;
    default:
        unreachable();
    }
}
