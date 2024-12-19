// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_TNG_FILE_HPP
#define CHEMFILES_TNG_FILE_HPP

#include <string>

#include <tng/tng_io_fwd.h>
#include <tng/tng_io.h>

#include "chemfiles/File.hpp"

namespace chemfiles {

/// Simple RAII capsule for `tng_trajectory_t`, handling the creation and
/// destruction of the file as needed.
class TNGFile final: public File {
public:
    TNGFile(std::string path, File::Mode mode);
    ~TNGFile() override;
    TNGFile(TNGFile&&) noexcept = default;
    TNGFile& operator=(TNGFile&&) = default;
    TNGFile(TNGFile const&) = delete;
    TNGFile& operator=(TNGFile const&) = delete;

    operator tng_trajectory_t() {return handle_;}

private:
    /// underlying pointer to the tng file
    tng_trajectory_t handle_;
};

/// Check a return code from a TNG function, and throw a `FileError` if the
/// status is not TNG_SUCCESS.
void check_tng_error(tng_function_status status, const std::string& function);

} // namespace chemfiles

#endif
