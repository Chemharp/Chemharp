// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>

#include "chemfiles/unreachable.hpp"
#include "chemfiles/error_fmt.hpp"

#include "chemfiles/File.hpp"
#include "chemfiles/files/PlainFile.hpp"

using namespace chemfiles;

#ifdef __CYGWIN__
    #include <sys/types.h>
    #define fseek64 fseek
    #define off64_t off_t
#elif defined(_MSC_VER)
    #define fseek64 _fseeki64
    #define off64_t __int64
#else
    // assume unix by default
    #include <sys/types.h>
    #define fseek64 fseeko
    #define off64_t off_t
    static_assert(_FILE_OFFSET_BITS == 64, "_FILE_OFFSET_BITS must be 64");
#endif

PlainFile::PlainFile(const std::string& path, File::Mode mode): TextFileImpl(path) {
    // We need to use binary mode when opening the file because we are storing
    // positions in the files relative to line ending positions. Using text
    // mode make the MSVC runtime convert lines ending and then all the values
    // return by tellpos are wrong.
    //
    // We can do this because we are dealing with line ending ourself.
    const char* openmode;

    switch (mode) {
    case File::READ:
        openmode = "rb";
        break;
    case File::APPEND:
        openmode = "a+b";
        break;
    case File::WRITE:
        openmode = "wb";
        break;
    default:
        unreachable();
    }

    file_ = std::fopen(path.c_str(), openmode);
    if (file_ == nullptr){
        throw file_error("could not open the file at '{}'", path);
    }
}

PlainFile::~PlainFile() {
    if (file_ != nullptr) {
        std::fclose(file_);
    }
}

void PlainFile::clear() noexcept {
    std::clearerr(file_);
}

void PlainFile::seek(uint64_t position) {
    static_assert(
        sizeof(uint64_t) == sizeof(off64_t),
        "uint64_t and off64_t do not have the same size"
    );
    auto status = fseek64(file_, static_cast<off64_t>(position), SEEK_SET);
    if (status != 0) {
        auto* message = std::strerror(errno);
        throw file_error("error while seeking file: {}", message);
    }
}

size_t PlainFile::read(char* data, size_t count) {
    auto result = std::fread(data, 1, count, file_);

    if (std::ferror(file_) != 0) {
        throw file_error("IO error while reading the file");
    }

    return result;
}

void PlainFile::write(const char* data, size_t count) {
    auto actual = std::fwrite(data, 1, count, file_);
    if (actual != count) {
        throw file_error("could not write data to the file at '{}'", this->path());
    }
}
