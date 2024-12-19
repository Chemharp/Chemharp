// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_FORMAT_TINKER_HPP
#define CHEMFILES_FORMAT_TINKER_HPP

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

/// Tinker XYZ file format.
///
/// The format is described here: http://chembytes.wikidot.com/tnk-tut00#toc2.
/// Additionally to the original Tinker XYZ format, chemfiles also supports the
/// files with unit cell data, as introduced in Tinker 6.3.
///
/// This format is associated with the .arc extension, but not the .xyz
/// extension, which is used for the standard XYZ format.
class TinkerFormat final: public TextFormat {
public:
    TinkerFormat(std::string path, File::Mode mode, File::Compression compression):
        TextFormat(std::move(path), mode, compression) {}
    TinkerFormat(std::shared_ptr<MemoryBuffer> memory, File::Mode mode, File::Compression compression) :
        TextFormat(std::move(memory), mode, compression) {}

    void read_next(Frame& frame) override;
    void write_next(const Frame& frame) override;
    optional<uint64_t> forward() override;
};

template<> const FormatMetadata& format_metadata<TinkerFormat>();

} // namespace chemfiles

#endif
