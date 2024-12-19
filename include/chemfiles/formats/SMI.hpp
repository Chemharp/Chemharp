// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_FORMAT_SMI_HPP
#define CHEMFILES_FORMAT_SMI_HPP

#include <cstddef>
#include <cstdint>

#include <map>
#include <stack>
#include <string_view>
#include <utility>
#include <vector>
#include <string>
#include <memory>
#include <unordered_map>

#include "chemfiles/File.hpp"
#include "chemfiles/Format.hpp"

#include "chemfiles/Residue.hpp"
#include "chemfiles/Connectivity.hpp"

#include "chemfiles/external/optional.hpp"

namespace chemfiles {
class Atom;
class Frame;
class Topology;
class MemoryBuffer;
class FormatMetadata;

/// SMI/OpenSMILES file format reader and writer.
class SMIFormat final: public TextFormat {
public:
    SMIFormat(std::string path, File::Mode mode, File::Compression compression):
        TextFormat(std::move(path), mode, compression) {}

    SMIFormat(std::shared_ptr<MemoryBuffer> memory, File::Mode mode, File::Compression compression) :
        TextFormat(std::move(memory), mode, compression) {}

    void read_next(Frame& frame) override;
    void write_next(const Frame& frame) override;
    optional<uint64_t> forward() override;

private:
    /// [for reading] adds an atom defined by `atom_name` to the topology
    Atom& add_atom(Topology& topology, std::string_view atom_name);

    /// [for reading] adds an atom defined by the string `smiles` starting at position i
    void process_property_list(Topology& topology, std::string_view smiles);

    /// [for reading] Opens and closes a ring with id `ring_id`
    void check_ring_(Topology& topology, size_t ring_id);

    /// [for reading] Stores location of a branching path
    std::stack<size_t, std::vector<size_t>> branch_point_;

    /// [for reading] Stores a mapping between a ring ID and the atom which
    /// starts the ring and a stored bond order
    std::unordered_map <size_t, std::pair<size_t, Bond::BondOrder>> rings_ids_;

    /// [for reading] The current atom being added (active atom)
    size_t current_atom_;

    /// [for reading] The previous active atom
    size_t previous_atom_;

    /// [for reading] The current bond order
    Bond::BondOrder current_bond_order_;

    /// [for reading] List of groups
    std::vector<Residue> residues_;

    /// [for reading] Should we connect the previous atom to the first atom
    /// [for writing] Should we add a '.' after the current molecule
    bool first_atom_;

    /// [for writing] stores the graph of the topology
    std::vector<std::vector<size_t>> adj_list_;

    /// [for writing] stores locations of ring closures. This is ordered
    /// to ensure rings get printed in numeric order (where possible)
    std::multimap<size_t, size_t> ring_stack_;

    /// [for writing] stores how many rings we need to close
    size_t ring_count_;

    /// [for writing] stores how many rings each atom is in
    std::unordered_map<size_t, size_t> ring_atoms_;
};

template<> const FormatMetadata& format_metadata<SMIFormat>();

} // namespace chemfiles

#endif
