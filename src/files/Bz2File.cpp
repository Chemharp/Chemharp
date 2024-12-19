// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cstdio>
#include <cstring>
#include <cstdint>
#include <cassert>
#include <limits>
#include <string>
#include <vector>
#include <functional>

#include <bzlib.h>

#include "chemfiles/File.hpp"
#include "chemfiles/error_fmt.hpp"

#include "chemfiles/files/MemoryBuffer.hpp"
#include "chemfiles/files/Bz2File.hpp"

using namespace chemfiles;

static unsigned checked_cast(uint64_t value) {
    if (value < std::numeric_limits<unsigned>::max()) {
        return static_cast<unsigned>(value);
    } else {
        throw file_error("{} is too big for unsigned in call to bzlib function", value);
    }
}

static void check(int status) {
    switch (status) {
    case BZ_OK:
    case BZ_RUN_OK:
    case BZ_FLUSH_OK:
    case BZ_FINISH_OK:
    case BZ_STREAM_END:
        return;
    case BZ_SEQUENCE_ERROR:
    case BZ_PARAM_ERROR:
        throw file_error("bzip2: bad call to bzlib (code: {}), this is a bug in chemfiles", status);
    case BZ_MEM_ERROR:
        throw file_error("bzip2: memory allocation failed (code: {})", status);
    case BZ_DATA_ERROR:
        throw file_error("bzip2: corrupted file (code: {})", status);
    case BZ_DATA_ERROR_MAGIC:
        throw file_error("bzip2: this file do not seems to be a bz2 file (code: {})", status);
    // These errors should not occur when using the stream API
    case BZ_CONFIG_ERROR:
        throw file_error("bzip2: mis-compiled bzlib (code: {})", status);
    case BZ_IO_ERROR:
    case BZ_UNEXPECTED_EOF:
    case BZ_OUTBUFF_FULL:
        throw file_error("bzip2: unexpected error from bzlib (code: {})", status);
    default:
        throw file_error("unknown error code from bzlib: {}", status);
    }
}

Bz2File::Bz2File(const std::string& path, File::Mode mode): TextFileImpl(path), mode_(mode), buffer_(8192) {
    std::memset(&stream_, 0, sizeof(bz_stream));

    const char* openmode = nullptr;
    if (mode == File::READ) {
        openmode = "rb";
        stream_end_ = BZ2_bzDecompressEnd;
        check(BZ2_bzDecompressInit(&stream_, 0, 0));
    } else if (mode == File::WRITE) {
        openmode = "wb";
        stream_end_ = BZ2_bzCompressEnd;
        check(BZ2_bzCompressInit(&stream_, 6, 0, 0));

        stream_.next_out = buffer_.data();
        stream_.avail_out = checked_cast(buffer_.size());
    } else if (mode == File::APPEND) {
        throw file_error("appending (open mode 'a') is not supported with bzip2 files");
    }

    file_ = std::fopen(path.c_str(), openmode);
    if (file_ == nullptr) {
        stream_end_(&stream_);
        throw file_error("could not open the file at '{}'", path);
    }
}

Bz2File::~Bz2File() {
    if (mode_ == File::WRITE) {
        try {
            compress_and_write(BZ_FINISH);
        } catch (...) {
            // not much we can do here
        }
    }

    stream_end_(&stream_);
    if (file_ != nullptr) {
        std::fclose(file_);
    }
}

size_t Bz2File::read(char* data, size_t count) {
    stream_.next_out = data;
    stream_.avail_out = checked_cast(count);

    while (stream_.avail_out != 0) {
        // read more compressed data from the file
        if (stream_.avail_in == 0 && (std::feof(file_) == 0)) {
            stream_.next_in = buffer_.data();
            stream_.avail_in = checked_cast(std::fread(buffer_.data(), 1, buffer_.size(), file_));

            if (std::ferror(file_) != 0) {
                throw file_error("IO error while reading bzip2 file");
            }
        }

        auto status = BZ2_bzDecompress(&stream_);

        if (status == BZ_STREAM_END) {
            return count - stream_.avail_out;
        } else {
            // Check for error
            check(status);
        }
    }
    return count;
}

void Bz2File::clear() noexcept {
    std::clearerr(file_);
}

void Bz2File::seek(uint64_t position) {
    assert(mode_ == File::READ);
    // Reset stream state
    stream_end_(&stream_);
    std::memset(&stream_, 0, sizeof(bz_stream));
    check(BZ2_bzDecompressInit(&stream_, 0, 0));

    // Dumb implementation, re-decompressing the file from the begining
    std::fseek(file_, 0, SEEK_SET);
    constexpr size_t BUFFSIZE = 4096;
    char buffer[BUFFSIZE];

    while (position > BUFFSIZE) {
        auto count = this->read(buffer, BUFFSIZE);
        assert(count == BUFFSIZE);
        position -= count;
    }

    auto count = this->read(buffer, static_cast<size_t>(position));
    assert(count == position);
    // silent "unused variable" when compiling without assertions
    (void)count;
}

void Bz2File::write(const char* data, size_t count) {
    stream_.next_in = const_cast<char*>(data);
    stream_.avail_in = checked_cast(count);
    compress_and_write(BZ_RUN);

    auto actual = count - stream_.avail_in;
    if (actual != count) {
        throw file_error("could not write data to the file at '{}'", this->path());
    }
}

void Bz2File::compress_and_write(int action) {
    int status = BZ_OK;
    do {
        status = BZ2_bzCompress(&stream_, action);

        if (stream_.avail_out == 0 || status == BZ_STREAM_END) {
            auto size = buffer_.size() - stream_.avail_out;
            auto written = std::fwrite(buffer_.data(), sizeof(uint8_t), size, file_);
            if (written != size) {
                throw file_error("error while writting data to bzip2 file");
            }

            stream_.next_out = buffer_.data();
            stream_.avail_out = checked_cast(buffer_.size());
        }

        check(status);

    } while (stream_.avail_in != 0 || (action == BZ_FINISH && status != BZ_STREAM_END));
}

// Get the full, potentially 64-bits large, value for total_out from the low and
// high 32-bits parts.
static uint64_t full_total_out(const bz_stream& stream) {
    uint64_t total_out_hi32 = static_cast<uint64_t>(stream.total_out_hi32) << 32;
    return total_out_hi32 + stream.total_out_lo32;
}

MemoryBuffer chemfiles::decompress_bz2(const char* src, size_t size) {
    // assume a 10% compression ratio, which should be plenty enough
    // (typical ratio is around 15-20%)
    auto output = MemoryBuffer(10 * size);

    bz_stream stream;
    stream.next_in = const_cast<char*>(src);
    stream.avail_in = checked_cast(size);
    stream.bzalloc = nullptr;
    stream.bzfree = nullptr;
    check(BZ2_bzDecompressInit(&stream, 0, 0));

    bool done = false;
    do {
        // if we need more space, resize the vector
        auto total_out = full_total_out(stream);
        if (total_out >= output.capacity()) {
            output.reserve_extra(output.capacity());
        }

	    stream.next_out = output.data_mut() + total_out;
        stream.avail_out = checked_cast(output.capacity() - total_out);

        auto status = BZ2_bzDecompress(&stream);
        if (status == BZ_STREAM_END) {
		    done = true;
        } else if (status != BZ_OK) {
		    BZ2_bzDecompressEnd(&stream);
            check(status);
	    }
    } while (!done);

    check(BZ2_bzDecompressEnd(&stream));

    auto total_out = full_total_out(stream);
    if (total_out >= output.capacity()) {
        // make sure the buffer always contains a terminal NULL
        output.reserve_extra(1);
    }
    output.set_size(checked_cast(total_out));
    return output;
}
