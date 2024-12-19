// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_FORMAT_MMTF_HPP
#define CHEMFILES_FORMAT_MMTF_HPP

#include <cstddef>
#include <cstdint>

#include <string>
#include <vector>
#include <memory>

#include <mmtf/structure_data.hpp>

#include "chemfiles/File.hpp"
#include "chemfiles/Format.hpp"

#include "chemfiles/UnitCell.hpp"
#include "chemfiles/external/span.hpp"

namespace chemfiles {
class Frame;
class Residue;
class Vector3D;
class MemoryBuffer;
class FormatMetadata;

/// MMTF file format reader and writer.
///
/// For multi-frame trajectories, we follow the PyMOL convention to use multiple
/// models for different frames
class MMTFFormat final: public Format {
public:
    MMTFFormat(std::string path, File::Mode mode, File::Compression compression);
    MMTFFormat(std::shared_ptr<MemoryBuffer> memory, File::Mode mode, File::Compression compression);

    ~MMTFFormat() override;

    void read_step(size_t step, Frame& frame) override;
    void read(Frame& frame) override;
    void write(const Frame& frame) override;
    size_t nsteps() override;

private:

    /// Perform the MMTF decoding steps
    void decode(const char* data, size_t size, const std::string& source);

    /// Add a model to a frame, increasing all the private indicies below
    void read_model(Frame& frame);

    /// Get the assembly that the current chain is a member of
    std::string find_assembly();

    /// Called by read_model to create new residues. Uses/updates groupIndex_
    Residue create_residue(const std::string& current_assembly, size_t group_type);

    /// Read a group from the MMTF structure, adding atoms and a residue to frame
    void read_group(Frame& frame, size_t group_type, Residue& residue, span<Vector3D> positions);

    /// Add inter residue bonds to the structure.
    void add_inter_residue_bonds(Frame& frame);

    /// Apply symmetry operations to the frame
    void apply_symmetry(Frame& frame);

    /// add a single residue to the structure_, using the data from the frame
    void add_residue_to_structure(const Frame& frame, const Residue& residue);

    /// A function to translate from the index in MMTF lists to an atom id
    /// suitable for chemfiles: starts at 0 for each model, and correspond to
    /// the initial atom index if it exists. This is used when reading.
    size_t atom_id(size_t mmtf_id);

    /// MMTF-CPP object holding all MacroMolecular information.
    /// This can be read and modified.
    mmtf::StructureData structure_;

    /// Location of MMTF file on disk. Only used if opened in write mode.
    std::string filename_;

    /// Current model being read. Ranges from [0, structure.numModels)
    size_t modelIndex_ = 0;

    /// Current chain being read. Ranges from [0, structure.numChains)
    size_t chainIndex_ = 0;

    /// Current group (residue) being read. Ranges from [0, structure.numGroups)
    size_t groupIndex_ = 0;

    /// Current atom being read. Ranges from [0, structure.numAtoms)
    size_t atomIndex_ = 0;

    /// Location of the inter-residue bonds being read.
    size_t interBondIndex_ = 0;

    /// Number of atoms before the current model.
    size_t atomSkip_ = 0;

    // Since MMTF uses model->chain->residue->atom as storage model, and
    // chemfiles do not enforce that residues contains contiguous atoms, the
    // atoms can be re-ordered when adding them to a MMTF structure. This vector
    // stores the association from chemfiles index to MMTF index to be able to
    // add the right bonds.
    //
    // This is only used when writing files, or reading files
    // with unordered or incomplete ids (e.g. reduced representation).
    std::vector<int32_t> new_atom_indexes_;

    /// Have we written the unitcell data?
    UnitCell unitcellForWrite_;
};

template<> const FormatMetadata& format_metadata<MMTFFormat>();

} // namespace chemfiles

#endif
