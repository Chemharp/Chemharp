// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cassert>
#include <cstdint>
#include <cstring>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <iterator>
#include <string_view>

#include <fmt/core.h>
#include <fmt/format.h>

#include "chemfiles/File.hpp"
#include "chemfiles/files/GzFile.hpp"
#include "chemfiles/files/XzFile.hpp"
#include "chemfiles/files/Bz2File.hpp"
#include "chemfiles/files/PlainFile.hpp"
#include "chemfiles/files/MemoryFile.hpp"
#include "chemfiles/files/MemoryBuffer.hpp"

#include "chemfiles/error_fmt.hpp"
#include "chemfiles/unreachable.hpp"

using namespace chemfiles;

TextFile::TextFile(std::string path, File::Mode mode, File::Compression compression):
    File(std::move(path), mode, compression),
    file_(nullptr),
    buffer_(8192, 0),
    line_start_(buffer_.data()),
    end_(buffer_.data() + buffer_.size())
{
    switch (compression) {
    case File::DEFAULT:
        file_ = std::make_unique<PlainFile>(this->path(), this->mode());
        break;
    case File::GZIP:
        file_ = std::make_unique<GzFile>(this->path(), this->mode());
        break;
    case File::BZIP2:
        file_ = std::make_unique<Bz2File>(this->path(), this->mode());
        break;
    case File::LZMA:
        file_ = std::make_unique<XzFile>(this->path(), this->mode());
        break;
    default:
        unreachable();
    }
}

TextFile::TextFile(std::shared_ptr<MemoryBuffer> memory, File::Mode mode, File::Compression compression):
    File("<in memory>", mode, File::Compression::DEFAULT),
    file_(nullptr),
    buffer_(8192, 0),
    line_start_(buffer_.data()),
    end_(buffer_.data() + buffer_.size())
{
    if (mode == File::APPEND) {
        throw file_error("cannot append (mode 'a') to a memory file");
    }

    if (compression != File::DEFAULT) {
        if (mode != File::READ) {
            throw file_error("writing to a compressed memory file is not supported");
        }

        memory->decompress(compression);
    }

    file_ = std::make_unique<MemoryFile>(std::move(memory), mode);
}

uint64_t TextFile::tellpos() const {
    assert(line_start_ >= buffer_.data());
    auto delta = buffer_initialized() ? static_cast<uint64_t>(line_start_ - buffer_.data()) : 0;
    return position_ + delta;
}

void TextFile::seekpos(uint64_t position) {
    got_impl_eof_ = false;
    eof_ = false;

    if (buffer_initialized()) {
        // use signed int64_t since the requested position can be smaller than
        // position_
        auto delta = static_cast<int64_t>(position) - static_cast<int64_t>(position_);
        if (0 <= delta && delta < static_cast<int64_t>(buffer_.size())) {
            // the new position is inside our buffer, no need to actually seek,
            // just reset the line_start_ to the corresponding position.
            line_start_ = buffer_.data() + delta;
            eof_ = false;
            return;
        }
    }

    // actually seek the file
    file_->seek(position);
    position_ = position;
    // mark buffer to be refilled
    buffer_[0] = '\0';
}

void TextFile::rewind() {
    clear();
    seekpos(0);
}

void TextFile::clear() {
    // Clear cached variables
    got_impl_eof_ = false;
    eof_ = false;
    // clear the file
    file_->clear();
}

bool TextFile::buffer_initialized() const {
    return buffer_[0] != '\0';
}

void TextFile::fill_buffer(size_t start) {
    auto count = buffer_.size() - start;
    if (buffer_initialized()) {
        position_ += count;
    }

    auto read_count = file_->read(buffer_.data() + start, count);
    if (read_count < count) {
        got_impl_eof_ = true;
        // Erase any remaining data in the buffer
        std::memset(buffer_.data() + start + read_count, 0, count - read_count);
    }
    line_start_ = buffer_.data();
}

std::string_view TextFile::readline() {
    // Initialize buffer if needed
    if (!buffer_initialized()) {
        fill_buffer(0);
    }

    if (eof_) {
        return "";
    }

    size_t length = 0;
    size_t windows_line = 0;
    while (true) {
        // How many characters are still in the buffer
        auto remainder = static_cast<size_t>(end_ - line_start_);
        // look for end of line character
        const auto *needle = std::memchr(line_start_ + length, '\n', remainder - length);
        const auto *newline = reinterpret_cast<const char*>(needle);

        if (newline != nullptr) {
            assert(line_start_ <= newline);
            length += static_cast<size_t>(newline - line_start_ + 1);
            // Check if we have a windows style line ending (\r\n)
            if (newline > buffer_.data() && newline[-1] == '\r') {
                windows_line = 1;
            }

            break;
        } else if (got_impl_eof_) {
            // no more data, we found the last line
            eof_ = true;

            // if line_start_ != end_ - 1, we still have a line in the buffer
            // (not terminated by a newline character)
            if (line_start_ != end_ - 1) {
                // We don't know the line length, but the buffer is terminated
                // with zeroes, so we can rely on std::strlen to create the
                // string_view.
                auto line = std::string_view(line_start_);
                line_start_ += line.length();
                return line;
            }
        }

        // no new line found in the current buffer, get more data
        if (remainder >= buffer_.size()) {
            // We did not found a new line in the whole buffer,
            // so we increase it size
            auto delta = line_start_ - buffer_.data();
            buffer_.resize(2 * buffer_.size(), 0);
            line_start_ = buffer_.data() + delta;
            end_ = buffer_.data() + buffer_.size();
        }

        // Move remaining data to the start of buffer
        std::memmove(buffer_.data(), line_start_, remainder);
        fill_buffer(remainder);
    }

    auto line = std::string_view(line_start_, length - windows_line - 1);
    line_start_ += length;

    return line;
}

void TextFile::vprint(fmt::string_view format, fmt::format_args args) {
    std::string buffer;
    buffer.reserve(128);
    fmt::vformat_to(std::back_inserter(buffer), format, args);
    if (buffer.size() == 0) {
        return;
    }
    file_->write(buffer.data(), buffer.size());
    position_ += buffer.size();
}

std::string TextFile::readall() {
    std::string buffer;
    buffer.resize(2048, '\0');
    size_t start = 0;
    while (true) {
        auto count = buffer.size() - start;
        auto read_count = file_->read(buffer.data() + start, count);
        start += read_count;
        if (read_count < count) {
            // Remove additional '\0' at the end
            buffer.resize(start);
            break;
        }
        buffer.resize(2 * buffer.size(), '\0');
    }

    return buffer;
}
