// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_FORMAT_GRO_HPP
#define CHEMFILES_FORMAT_GRO_HPP

#include <cstdint>
#include <map>
#include <string>
#include <memory>
#include <utility>

#include "chemfiles/File.hpp"
#include "chemfiles/Format.hpp"

#include "chemfiles/Residue.hpp"
#include "chemfiles/external/optional.hpp"

namespace chemfiles {
class Frame;
class MemoryBuffer;
class FormatMetadata;

/// GRO file format reader and writer.
class GROFormat final: public TextFormat {
public:
    GROFormat(std::string path, File::Mode mode, File::Compression compression):
        TextFormat(std::move(path), mode, compression) {}

    GROFormat(std::shared_ptr<MemoryBuffer> memory, File::Mode mode, File::Compression compression) :
        TextFormat(std::move(memory), mode, compression) {}

    void read_next(Frame& frame) override;
    void write_next(const Frame& frame) override;
    optional<uint64_t> forward() override;

private:
    /// Map of residues, indexed by residue id.
    std::map<int64_t, Residue> residues_;
};

template<> const FormatMetadata& format_metadata<GROFormat>();

} // namespace chemfiles

#endif
