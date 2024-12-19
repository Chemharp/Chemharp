// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_FORMAT_SDF_HPP
#define CHEMFILES_FORMAT_SDF_HPP

#include <cstdint>
#include <string>
#include <memory>
#include <utility>

#include "chemfiles/File.hpp"
#include "chemfiles/Format.hpp"
#include "chemfiles/external/optional.hpp"

namespace chemfiles {
class Frame;
class MemoryBuffer;
class FormatMetadata;

/// SDF file format reader and writer.
class SDFFormat final: public TextFormat {
public:
    SDFFormat(std::string path, File::Mode mode, File::Compression compression):
        TextFormat(std::move(path), mode, compression) {}

    SDFFormat(std::shared_ptr<MemoryBuffer> memory, File::Mode mode, File::Compression compression) :
        TextFormat(std::move(memory), mode, compression) {}

    void read_next(Frame& frame) override;
    void write_next(const Frame& frame) override;
    optional<uint64_t> forward() override;
};

template<> const FormatMetadata& format_metadata<SDFFormat>();

} // namespace chemfiles

#endif
