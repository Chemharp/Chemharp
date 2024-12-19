// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_MEMORY_FILES_HPP
#define CHEMFILES_MEMORY_FILES_HPP

#include <cstdint>
#include <cstddef>
#include <memory>
#include <utility>

#include "chemfiles/File.hpp"

namespace chemfiles {
class MemoryBuffer;

/// Simple TextFile implementation, that wraps memory for reading only.
class MemoryFile final: public TextFileImpl {
public:
    /// Open `memory` as though it were a file in mode `mode`. No copy of `memory` is
    /// made and the original object **MUST** not be freed until this object is destroyed
    MemoryFile(std::shared_ptr<MemoryBuffer> memory, File::Mode mode)
        : TextFileImpl(""), buffer_(std::move(memory)), mode_(mode)
    {}

    size_t read(char* data, size_t count) override;
    void write(const char* data, size_t count) override;

    void clear() noexcept override {}
    void seek(uint64_t position) override;

private:
    /// Current reading location
    size_t current_location_ = 0;

    /// An input-out stream that performs all the read operations, shared with
    /// the trajectory
    std::shared_ptr<MemoryBuffer> buffer_;

    /// Is this for reading or writing?
    File::Mode mode_;
};

} // namespace chemfiles

#endif
