// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_FORMAT_XYZ_HPP
#define CHEMFILES_FORMAT_XYZ_HPP

#include <cstddef>
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

/// XYZ file format reader and writer.
///
/// This class also support the extended XYZ specification, as defined in
/// [ASE](https://wiki.fysik.dtu.dk/ase/ase/io/formatoptions.html#extxyz)
class XYZFormat final: public TextFormat {
public:
    XYZFormat(std::string path, File::Mode mode, File::Compression compression):
        TextFormat(std::move(path), mode, compression){}

    XYZFormat(std::shared_ptr<MemoryBuffer> memory, File::Mode mode, File::Compression compression) :
        TextFormat(std::move(memory), mode, compression){}

    void read_next(Frame& frame) override;
    void write_next(const Frame& frame) override;
    optional<uint64_t> forward() override;

private:
    // used to give better error message in `forward`, this refers to the
    // current step being checked.
    size_t current_forward_step_ = 0;
};

template<> const FormatMetadata& format_metadata<XYZFormat>();

} // namespace chemfiles

#endif
