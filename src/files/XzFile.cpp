// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cassert>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>
#include <limits>

#include <lzma.h>

#include "chemfiles/File.hpp"
#include "chemfiles/error_fmt.hpp"

#include "chemfiles/files/XzFile.hpp"
#include "chemfiles/files/MemoryBuffer.hpp"

using namespace chemfiles;

static size_t checked_cast(uint64_t value) {
    if (value < std::numeric_limits<size_t>::max()) {
        return static_cast<size_t>(value);
    } else {
        throw file_error("{} is too big for size_t in call to xz function", value);
    }
}

static void check(lzma_ret code) {
    int int_code = static_cast<int>(code);
    switch (code) {
    case LZMA_OK:
    case LZMA_STREAM_END:
        return;
    case LZMA_MEM_ERROR:
    case LZMA_MEMLIMIT_ERROR:
        throw file_error("lzma: memory allocation failed (code: {})", int_code);
    case LZMA_FORMAT_ERROR:
        throw file_error("lzma: input not in .xz format (code: {})", int_code);
    case LZMA_OPTIONS_ERROR:
        throw file_error(
            "lzma: unsupported compression options (code: {})", int_code);
    case LZMA_DATA_ERROR:
        throw file_error("lzma: compressed file is corrupted (code: {})", int_code);
    case LZMA_BUF_ERROR:
        throw file_error(
            "lzma: compressed file is truncated or corrupted (code: {})", int_code);
    case LZMA_UNSUPPORTED_CHECK:
        throw file_error("lzma: specified integrity check is not supported (code: {})", int_code);
    case LZMA_PROG_ERROR:
        throw file_error("lzma: programming error (this is a bug) (code: {})", int_code);
    default:
        throw file_error("lzma: unknown error (code: {})", int_code);
    }
}

static void open_stream_read(lzma_stream* stream) {
    auto flags = LZMA_TELL_UNSUPPORTED_CHECK | LZMA_CONCATENATED;
    auto memory_limit = std::numeric_limits<uint64_t>::max();
    check(lzma_stream_decoder(stream, memory_limit, flags));
}

XzFile::XzFile(const std::string& path, File::Mode mode): TextFileImpl(path), mode_(mode), buffer_(8192) {
    const char* openmode = nullptr;
    if (mode == File::READ) {
        openmode = "rb";
        open_stream_read(&stream_);
    } else if (mode == File::WRITE) {
        openmode = "wb";
        check(lzma_easy_encoder(&stream_, 6, LZMA_CHECK_CRC64));

        stream_.next_out = buffer_.data();
        stream_.avail_out = buffer_.size();
    } else if (mode == File::APPEND) {
        throw file_error("appending (open mode 'a') is not supported with xz files");
    }

    file_ = std::fopen(path.c_str(), openmode);
    if (file_ == nullptr) {
        lzma_end(&stream_);
        throw file_error("could not open the file at '{}'", path);
    }
}

XzFile::~XzFile() {
    if (mode_ == File::WRITE) {
        try {
            compress_and_write(LZMA_FINISH);
        } catch (...) {
            // not much we can do here ...
        }
    }

    lzma_end(&stream_);
    if (file_ != nullptr) {
        std::fclose(file_);
    }
}

size_t XzFile::read(char* data, size_t count) {
    auto action = LZMA_RUN;

    stream_.next_out = reinterpret_cast<uint8_t*>(data);
    stream_.avail_out = count;

    while (stream_.avail_out != 0) {
        // read more compressed data from the file
        if (stream_.avail_in == 0 && (std::feof(file_) == 0)) {
            stream_.next_in = buffer_.data();
            stream_.avail_in = std::fread(buffer_.data(), 1, buffer_.size(), file_);

            if (std::ferror(file_) != 0) {
                throw file_error("IO error while reading xz file");
            }
        }

        if (std::feof(file_) != 0) {
            action = LZMA_FINISH;
        }

        auto status = lzma_code(&stream_, action);

        if (status == LZMA_STREAM_END) {
            return count - stream_.avail_out;
        } else {
            // Check for error
            check(status);
        }
    }
    return count;
}

void XzFile::clear() noexcept {
    std::clearerr(file_);
}

void XzFile::seek(uint64_t position) {
    assert(mode_ == File::READ);
    // Reset stream state
    lzma_end(&stream_);
    stream_ = LZMA_STREAM_INIT;
    open_stream_read(&stream_);

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

void XzFile::write(const char* data, size_t count) {
    stream_.next_in = reinterpret_cast<const uint8_t*>(data);
    stream_.avail_in = count;
    compress_and_write(LZMA_RUN);

    auto actual = count - stream_.avail_in;
    if (actual != count) {
        throw file_error("could not write data to the file at '{}'", this->path());
    }
}

void XzFile::compress_and_write(lzma_action action) {
    lzma_ret status = LZMA_OK;
    do {
        status = lzma_code(&stream_, action);

        if (stream_.avail_out == 0 || status == LZMA_STREAM_END) {
            auto size = buffer_.size() - stream_.avail_out;
            auto written = std::fwrite(buffer_.data(), sizeof(uint8_t), size, file_);
            if (written != size) {
                throw file_error("error while writting data to xz file");
            }

            stream_.next_out = buffer_.data();
            stream_.avail_out = buffer_.size();
        }

        check(status);

    } while (stream_.avail_in != 0 || (action == LZMA_FINISH && status != LZMA_STREAM_END));
}

MemoryBuffer chemfiles::decompress_xz(const char* src, size_t size) {
    // assume a 10% compression ratio, which should be plenty enough
    // (typical ratio is around 15-20%)
    auto output = MemoryBuffer(10 * size);

    lzma_stream stream = LZMA_STREAM_INIT;
    stream.next_in = reinterpret_cast<const uint8_t*>(src);
    stream.avail_in = size;
    stream.total_out = 0;

    open_stream_read(&stream);

    bool done = false;
    do {
	    // if we need more space, resize the vector
	    if (stream.total_out >= output.capacity()) {
		    output.reserve_extra(output.capacity());
	    }

	    stream.next_out = reinterpret_cast<uint8_t*>(output.data_mut() + stream.total_out);
    	stream.avail_out = output.capacity() - checked_cast(stream.total_out);

        auto status = lzma_code(&stream, LZMA_FINISH);
        if (status == LZMA_STREAM_END) {
		    done = true;
        } else if (status != LZMA_OK) {
		    lzma_end(&stream);
            check(status);
	    }
    } while (!done);

    lzma_end(&stream);

    if (stream.total_out >= output.capacity()) {
        // make sure the buffer always contains a terminal NULL
        output.reserve_extra(1);
    }
    output.set_size(checked_cast(stream.total_out));
    return output;
}
