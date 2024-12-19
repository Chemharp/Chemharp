// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_FORMAT_MOL2_HPP
#define CHEMFILES_FORMAT_MOL2_HPP

#include <cstddef>
#include <cstdint>

#include <string>
#include <memory>
#include <unordered_map>
#include <utility>

#include "chemfiles/File.hpp"
#include "chemfiles/Format.hpp"
#include "chemfiles/external/optional.hpp"

namespace chemfiles {
class Frame;
class Residue;
class MemoryBuffer;
class FormatMetadata;

/// Tripos MOL2 format reader and writer.
class MOL2Format final: public TextFormat {
public:
    MOL2Format(std::string path, File::Mode mode, File::Compression compression):
        TextFormat(std::move(path), mode, compression) {}

    MOL2Format(std::shared_ptr<MemoryBuffer> memory, File::Mode mode, File::Compression compression) :
        TextFormat(std::move(memory), mode, compression) {}

    void read_next(Frame& frame) override;
    void write_next(const Frame& frame) override;
    optional<uint64_t> forward() override;

private:
    // Read Atoms
    void read_atoms(Frame& frame, size_t natoms, bool charges);

    // Read Bonds
    void read_bonds(Frame& frame, size_t nbonds);

    /// Map of residues, indexed by residue id.
    std::unordered_map<int64_t, Residue> residues_;
};

template<> const FormatMetadata& format_metadata<MOL2Format>();

} // namespace chemfiles

#endif
