// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cstddef>
#include <cstdint>

#include <memory>
#include <algorithm>

#include "chemfiles/File.hpp"
#include "chemfiles/files/MemoryFile.hpp"
#include "chemfiles/files/MemoryBuffer.hpp"  // IWYU pragma: keep

#include "chemfiles/error_fmt.hpp"

using namespace chemfiles;

void MemoryFile::seek(uint64_t position) {
    if (mode_ != File::READ) {
        throw file_error("cannot seek a memory file unless it is opened in read mode");
    }

    if (position > buffer_->size()) {
        current_location_ =  buffer_->size();
    }

    current_location_ = static_cast<size_t>(position);
}

size_t MemoryFile::read(char* data, size_t count) {
    if (mode_ != File::READ) {
        throw file_error("cannot read a memory file unless it is opened in read mode");
    }

    if (current_location_ >= buffer_->size()) {
        return 0;
    }

    auto amount_to_read = count + current_location_ <= buffer_->size() ?
        count : buffer_->size() - current_location_;

    const char* start = buffer_->data() + current_location_;
    std::copy(start, start + amount_to_read, data);
    current_location_ += amount_to_read;

    return amount_to_read;
}

void MemoryFile::write(const char* data, size_t count) {
    if (mode_ != File::WRITE) {
        throw file_error("cannot write to a memory file unless it is opened in write mode");
    }

    buffer_->write(data, count);
}
