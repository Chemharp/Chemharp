// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_MEMORY_FILES_HPP
#define CHEMFILES_MEMORY_FILES_HPP

#include "chemfiles/File.hpp"
#include "chemfiles/files/MemoryBuffer.hpp"

#include <istream> // Needed for clang
#include <ostream> // Needed for clang
#include <iosfwd>

namespace chemfiles {

/// Simple TextFile implementation, that wraps memory for reading only.
class MemoryFile final: public TextFileImpl {

public:
    /// Open `memory` as though it were a file in mode `mode`. No copy of `memory` is
    /// made and the original object **MUST** not be freed until this object is destroyed
    MemoryFile(MemoryBuffer& memory, File::Mode mode)
        : TextFileImpl(""), current_location_(0), data_(memory), mode_(mode)
    {}

    size_t read(char* data, size_t count) override;
    void write(const char* data, size_t count) override;

    void clear() noexcept override;
    void seek(uint64_t position) override;

private:

    /// Current reading location
    size_t current_location_;

    /// An input-out stream that performs all the read operations
    MemoryBuffer& data_;

    /// Is this for reading or writing?
    File::Mode mode_;
};

} // namespace chemfiles

#endif
