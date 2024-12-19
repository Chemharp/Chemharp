// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <string>
#include <utility>

#include <tng/tng_io.h>

#include "chemfiles/File.hpp"
#include "chemfiles/utils.hpp"
#include "chemfiles/error_fmt.hpp"

#include "chemfiles/files/TNGFile.hpp"

using namespace chemfiles;

#define STRING_0(x) #x
#define STRING(x) STRING_0(x)
#define CHECK(x) check_tng_error((x), (STRING(x)))

TNGFile::TNGFile(std::string path, File::Mode mode): File(std::move(path), mode, File::DEFAULT), handle_(nullptr) {
    CHECK(tng_util_trajectory_open(this->path().c_str(), mode, &handle_));

    if (mode == File::READ) {
        auto status = tng_file_headers_read(handle_, TNG_USE_HASH);
        if (status != TNG_SUCCESS) {
            tng_util_trajectory_close(&handle_);
            throw file_error("could not open the file at '{}'", this->path());
        }
    }

    if (mode == File::WRITE || mode == File::APPEND) {
        CHECK(tng_last_program_name_set(handle_, "chemfiles"));
        auto user = user_name();
        CHECK(tng_last_user_name_set(handle_, user.c_str()));
        auto host = hostname();
        CHECK(tng_last_computer_name_set(handle_, host.c_str()));
        // TODO: last pgp signature?

        if (mode == File::WRITE) {
            CHECK(tng_first_program_name_set(handle_, "chemfiles"));
            CHECK(tng_first_user_name_set(handle_, user.c_str()));
            CHECK(tng_first_computer_name_set(handle_, host.c_str()));
            // TODO: first pgp signature?
        }

        CHECK(tng_file_headers_write(handle_, TNG_USE_HASH));
    }
}

TNGFile::~TNGFile() {
    tng_util_trajectory_close(&handle_);
}

void chemfiles::check_tng_error(tng_function_status status, const std::string& function) {
    switch (status) {
    case TNG_FAILURE:
        throw file_error("error while calling {} in the TNG library", function);
    case TNG_CRITICAL:
        throw file_error("critical error while calling {} in the TNG library", function);
    case TNG_SUCCESS:
        // Do nothing, this is good
        break;
    default:
        throw file_error("unknown status code from TNG library: {}", static_cast<int>(status));
    }
}
