// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_FORMAT_PDB_HPP
#define CHEMFILES_FORMAT_PDB_HPP

#include <cstddef>
#include <cstdint>

#include <map>
#include <string>
#include <string_view>
#include <utility>
#include <vector>
#include <memory>

#include "chemfiles/File.hpp"
#include "chemfiles/Format.hpp"

#include "chemfiles/Residue.hpp"
#include "chemfiles/external/optional.hpp"

namespace chemfiles {
class Frame;
class MemoryBuffer;
class FormatMetadata;

/// Full identification of residues, including everything needed
struct FullResidueId {
    /// Chain identifier
    char chain;
    /// Residue id
    int64_t resid;
    /// Residue name
    std::string resname;
    /// Insertion code of the residue
    char insertion_code;
};

bool operator==(const FullResidueId& lhs, const FullResidueId& rhs);
bool operator<(const FullResidueId& lhs, const FullResidueId& rhs);

/// PDB file format reader and writer.
///
/// For multi-frame trajectories, we support both the convention from VMD to
/// use multiple `END` records separating the steps; or the use of multiple
/// `MODEL`/`ENDMODEL` pairs.
class PDBFormat final: public TextFormat {
public:
    PDBFormat(std::string path, File::Mode mode, File::Compression compression):
        TextFormat(std::move(path), mode, compression) {}

    PDBFormat(std::shared_ptr<MemoryBuffer> memory, File::Mode mode, File::Compression compression) :
        TextFormat(std::move(memory), mode, compression) {}

    ~PDBFormat() override;

    void read_next(Frame& frame) override;
    void write_next(const Frame& frame) override;
    optional<uint64_t> forward() override;

    // Connect residues based on a predefined table
    static void link_standard_residue_bonds(Frame& frame);

private:
    // Read HEADER record
    void read_HEADER(Frame& frame, std::string_view line);
    // Read TITLE record
    void read_TITLE(Frame& frame, std::string_view line);
    // Read CRYST1 record
    void read_CRYST1(Frame& frame, std::string_view line);
    // Read ATOM and HETATM records
    void read_ATOM(Frame& frame, std::string_view line, bool is_hetatm);
    // Read secondary structure records. All push secinfo_ vector if line is valid
    void read_HELIX(std::string_view line);
    // reads SHEET and TURN records. i_start and i_end are the indicies of the chain ids
    void read_secondary(std::string_view line, size_t i_start, size_t i_end, std::string_view record);
    // Read CONECT record
    void read_CONECT(Frame& frame, std::string_view line);
    // Runs when a chain is terminated to update residue information
    void chain_ended(Frame& frame);

    /// Residue information in the current step
    std::map<FullResidueId, Residue> residues_;
    /// Number of models written/read to the file.
    size_t models_ = 0;
    /// List of all atom offsets. This maybe pushed in read_ATOM or if a TER
    /// record is found. It is reset every time a frame is read.
    std::vector<size_t> atom_offsets_;
    /// Did we wrote a frame to the file? This is used to check whether we need
    /// to write a final `END` record in the destructor
    bool written_ = false;
    /// Store secondary structure information. Keys are the
    /// starting residue of the secondary structure, and values are pairs
    /// containing the ending residue and a string which is a written
    /// description of the secondary structure
    std::map<FullResidueId, std::pair<FullResidueId, std::string>> secinfo_;
    /// This will be nullopt when no secondary structure information should be
    /// read. Else It is set to the final residue of a secondary structure and
    /// the text description which should be set.
    optional<std::pair<FullResidueId, std::string>> current_secinfo_;
};

template<> const FormatMetadata& format_metadata<PDBFormat>();

} // namespace chemfiles

#endif
