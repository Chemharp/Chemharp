// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>

#include "chemfiles/File.hpp"
#include "chemfiles/files/PlainFile.hpp"
#include "chemfiles/string_view.hpp"
#include "chemfiles/error_fmt.hpp"

using namespace chemfiles;

#ifdef __CYGWIN__
    #define fseek64 fseek
#elif defined(_MSC_VER)
    #define fseek64 _fseeki64
#else
    // assume unix by default
    #define fseek64 fseeko
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
// TODO: update to a more recent emscripten, with support for 64-bit fseek
#if !defined(__EMSCRIPTEN__)
    static_assert(
        sizeof(uint64_t) == sizeof(off_t),
        "uint64_t and off_t do not have the same size"
    );
#endif
    auto status = fseek64(file_, static_cast<off_t>(position), SEEK_SET);
    if (status != 0) {
        auto message = std::strerror(errno);
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
