// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_FORMAT_MMCIF_HPP
#define CHEMFILES_FORMAT_MMCIF_HPP

#include <cstddef>
#include <cstdint>

#include <map>
#include <string>
#include <utility>
#include <vector>
#include <memory>

#include "chemfiles/File.hpp"
#include "chemfiles/Format.hpp"

#include "chemfiles/Residue.hpp"
#include "chemfiles/UnitCell.hpp"

namespace chemfiles {
class Frame;
class MemoryBuffer;
class FormatMetadata;

/// mmCIF Crystallographic Information Framework for MacroMolecules
/// reader and writer.
class mmCIFFormat final: public Format {
public:
    mmCIFFormat(std::string path, File::Mode mode, File::Compression compression) :
        file_(std::move(path), mode, compression), models_(0), atoms_(0) {
        init_();
    }

    mmCIFFormat(std::shared_ptr<MemoryBuffer> memory, File::Mode mode, File::Compression compression) :
        file_(std::move(memory), mode, compression), models_(0), atoms_(0) {
        init_();
    }

    void read_step(size_t step, Frame& frame) override;
    void read(Frame& frame) override;
    void write(const Frame& frame) override;
    size_t nsteps() override;
private:
    /// Initialize important variables
    void init_();
    /// Underlying file representation
    TextFile file_;
    /// Map of STAR records to their index
    std::map<std::string, size_t> atom_site_map_;
    /// Vector with all the residues.
    std::vector<Residue> residues_;
    /// Map of residue indexes, indexed by residue id and chainid. We use an indirection to keep the residue order (and don't sort them with the map id).
    std::map<std::pair<std::string, int64_t>, size_t> map_residues_indexes;
    /// Storing the positions of all the steps in the file, so that we can
    /// just `seekpos` them instead of reading the whole step.
    std::vector<uint64_t> steps_positions_;
    /// The cell for all frames
    UnitCell cell_;
    /// Number of models written to the file.
    size_t models_;
    /// Number of atoms written to the file.
    size_t atoms_;
    /// Frame properties need to be stored
    std::string name_;
    /// The PDB icode, if any
    std::string pdb_idcode_;
};

template<> const FormatMetadata& format_metadata<mmCIFFormat>();

} // namespace chemfiles

#endif
