// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license
#include <cstdio>
#include <cstdint>
#include <limits>
#include <string>

#define ZLIB_CONST
#include <zconf.h>
#include <zlib.h>

#include "chemfiles/File.hpp"
#include "chemfiles/error_fmt.hpp"
#include "chemfiles/unreachable.hpp"

#include "chemfiles/files/MemoryBuffer.hpp"
#include "chemfiles/files/GzFile.hpp"

using namespace chemfiles;

static unsigned checked_cast(size_t value) {
    if (value < std::numeric_limits<unsigned>::max()) {
        return static_cast<unsigned>(value);
    } else {
        throw file_error("{} is too big for unsigned in call to zlib function", value);
    }
}

GzFile::GzFile(const std::string& path, File::Mode mode): TextFileImpl(path) {
    const char* openmode;
    switch (mode) {
    case File::READ:
        openmode = "rb";
        break;
    case File::WRITE:
        openmode = "wb7";
        break;
    case File::APPEND:
        openmode = "ab7";
        break;
    default:
        unreachable();
    }

    file_ = gzopen64(path.c_str(), openmode);
    if (file_ == nullptr) {
        throw file_error("could not open the file at '{}'", path);
    }
}

GzFile::~GzFile() {
    if (file_ != nullptr) {
        gzclose(file_);
    }
}

size_t GzFile::read(char* data, size_t count) {
    auto read = gzread(file_, data, checked_cast(count));
    const auto* error = check_error();
    if (read == -1 || error != nullptr) {
        throw file_error("error while reading gziped file: {}", error);
    }
    return static_cast<size_t>(read);
}

void GzFile::write(const char* data, size_t count) {
    auto actual = gzwrite(file_, data, checked_cast(count));
    const auto* error = check_error();
    if (actual == 0 || error != nullptr) {
        throw file_error("error while writting to gziped file: {}", error);
    }
    if (static_cast<size_t>(actual) != count) {
        throw file_error("could not write data to the file at '{}'", this->path());
    }
}

const char* GzFile::check_error() const {
    int status = Z_OK;
    const auto* message = gzerror(file_, &status);
    return status == Z_OK ? nullptr : message;
}

void GzFile::clear() noexcept {
    gzclearerr(file_);
}

void GzFile::seek(uint64_t position) {
    static_assert(
        sizeof(uint64_t) == sizeof(z_off64_t),
        "uint64_t and z_off64_t do not have the same size"
    );
    auto status = gzseek64(file_, static_cast<z_off64_t>(position), SEEK_SET);
    if (status == -1) {
        const auto* message = check_error();
        throw file_error("error while seeking gziped file: {}", message);
    }
}

MemoryBuffer chemfiles::decompress_gz(const char* src, size_t size) {
    // assume a 10% compression ratio, which should be plenty enough
    // (typical ratio is around 15-20%)
    auto output = MemoryBuffer(10 * size);

    z_stream stream;
    stream.next_in = reinterpret_cast<const Bytef*>(src);
    stream.avail_in = checked_cast(size);
    stream.total_out = 0;
    stream.zalloc = nullptr;
    stream.zfree = nullptr;

    // the second parameter is set to 15 (use the largest window possible) + 32
    // (detect header and check between gzip or zlib header)
    auto status = inflateInit2(&stream, 15 + 32);
    if (status != Z_OK) {
	    throw file_error("error creating gz stream: {}", stream.msg);
    }

    bool done = false;
    do {
        // if we need more space, resize the vector
        if (stream.total_out >= output.capacity()) {
            output.reserve_extra(output.capacity());
        }

	    stream.next_out = reinterpret_cast<Bytef*>(output.data_mut() + stream.total_out);
        stream.avail_out = checked_cast(output.capacity() - stream.total_out);

        status = inflate(&stream, Z_SYNC_FLUSH);
        if (status == Z_STREAM_END) {
		    done = true;
        } else if (status != Z_OK) {
		    inflateEnd(&stream);
            throw file_error("error inflating gziped memory: {}", stream.msg);
	    }
    } while (!done);

    status = inflateEnd(&stream);
    if (status != Z_OK) {
	    throw file_error("error finishing gz stream: {}", stream.msg);
    }

    if (stream.total_out >= output.capacity()) {
        // make sure the buffer always contains a terminal NULL
        output.reserve_extra(1);
    }
    output.set_size(stream.total_out);
    return output;
}
