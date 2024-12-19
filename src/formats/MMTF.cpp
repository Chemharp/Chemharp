// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cstddef>
#include <cstdint>
#include <cassert>

#include <algorithm>
#include <array>
#include <iterator>
#include <string>
#include <utility>
#include <vector>
#include <memory>
#include <exception>
#include <unordered_set>

#include <mmtf/errors.hpp>
#include <mmtf/structure_data.hpp>
#include <mmtf/decoder.hpp>
#include <mmtf/encoder.hpp>
#include <mmtf/export_helpers.hpp>

#include "chemfiles/types.hpp"
#include "chemfiles/warnings.hpp"
#include "chemfiles/error_fmt.hpp"
#include "chemfiles/unreachable.hpp"
#include "chemfiles/external/span.hpp"
#include "chemfiles/external/optional.hpp"

#include "chemfiles/File.hpp"
#include "chemfiles/Atom.hpp"
#include "chemfiles/Frame.hpp"
#include "chemfiles/Residue.hpp"
#include "chemfiles/UnitCell.hpp"
#include "chemfiles/Topology.hpp"
#include "chemfiles/Property.hpp"
#include "chemfiles/Connectivity.hpp"
#include "chemfiles/FormatMetadata.hpp"

#include "chemfiles/files/MemoryBuffer.hpp"

#include "chemfiles/formats/MMTF.hpp"

using namespace chemfiles;

static int8_t bond_order_to_mmtf(Bond::BondOrder order);
static Bond::BondOrder bond_order_to_chemfiles(int32_t order);
static void set_secondary(Residue& residue, int32_t code);

template<> const FormatMetadata& chemfiles::format_metadata<MMTFFormat>() {
    static FormatMetadata metadata;
    metadata.name = "MMTF";
    metadata.extension = ".mmtf";
    metadata.description = "MMTF (RCSB Protein Data Bank) binary format";
    metadata.reference = "https://mmtf.rcsb.org/";

    metadata.read = true;
    metadata.write = true;
    metadata.memory = true;

    metadata.positions = true;
    metadata.velocities = false;
    metadata.unit_cell = true;
    metadata.atoms = true;
    metadata.bonds = true;
    metadata.residues = true;
    return metadata;
}

MMTFFormat::MMTFFormat(std::string path, File::Mode mode, File::Compression compression) {
    if (mode == File::READ) {
        auto file = TextFile(std::move(path), mode, compression);
        auto buffer = file.readall();
        decode(buffer.data(), buffer.size(), file.path());
        if (!mmtf::isDefaultValue(structure_.atomIdList)) {
            // If ids are not ordered or are missing consecutive values the atoms
            // have to be re-ordered.
            bool isValidIdOrder = false;
            if (structure_.atomIdList[0] == 1) {
                // hack `is_sorted` to check that all ids are consecutive and sorted
                isValidIdOrder = std::is_sorted(
                    new_atom_indexes_.begin(), new_atom_indexes_.end(),
                    [](const int32_t& lhs, const int32_t& rhs) { return (lhs + 1) == rhs; });
            }
            if (!isValidIdOrder) {
                new_atom_indexes_ = structure_.atomIdList;
                std::sort(new_atom_indexes_.begin(), new_atom_indexes_.end());
            }
        }
    } else if (mode == File::WRITE) {
        filename_ = std::move(path); // We really don't need to do anything, yet
    } else if (mode == File::APPEND) {
        throw file_error("append mode ('a') is not supported for the MMTF format");
    }
}

MMTFFormat::MMTFFormat(std::shared_ptr<MemoryBuffer> memory, File::Mode mode, File::Compression compression) {
    if (mode == File::WRITE) {
        throw format_error("the MMTF format cannot write to memory");
    }

    memory->decompress(compression);
    decode(memory->data(), memory->size(), "memory");
}

void MMTFFormat::decode(const char* data, size_t size, const std::string& source) {
    try {
        mmtf::decodeFromBuffer(structure_, data, size);
    } catch (const mmtf::DecodeError& e) { // rethrow as a chemfiles error
        throw format_error("error while decoding MMTF from {}: '{}'", source, e.what());
    }

    if (!structure_.hasConsistentData()) {
        throw format_error("issue with data from '{}', please ensure it is valid MMTF file", source);
    }
}

MMTFFormat::~MMTFFormat() {
    if (!filename_.empty()) {
        try {
            mmtf::compressGroupList(structure_);
            encodeToFile(structure_, filename_);
        } catch (const std::exception& e) {
            warning("MMTF writer", "error while finishing writing to {}: {}", filename_, e.what());
        } catch (...) {
            // ignore exceptions in destructor
        }
    }
}

size_t MMTFFormat::nsteps() {
    return static_cast<size_t>(structure_.numModels);
}

void MMTFFormat::read_step(const size_t step, Frame& frame) {
    modelIndex_ = 0;
    chainIndex_ = 0;
    groupIndex_ = 0;
    atomIndex_  = 0;
    atomSkip_   = 0;
    interBondIndex_ = 0;

    // Fast-forward, keeping all indexes updated
    while(modelIndex_ != step) {
        auto chainsPerModel = static_cast<size_t>(structure_.chainsPerModel[modelIndex_]);
        for (size_t j = 0; j < chainsPerModel; ++j) {
            auto groupsPerChain = static_cast<size_t>(structure_.groupsPerChain[chainIndex_]);
            for (size_t k = 0; k < groupsPerChain; ++k) {
                auto groupType = static_cast<size_t>(structure_.groupTypeList[groupIndex_]);
                const auto& group = structure_.groupList[groupType];
                auto atomCount = group.atomNameList.size();
                atomIndex_ += atomCount;
                groupIndex_++;
            }
            chainIndex_++;
        }
        modelIndex_++;
    }

    atomSkip_ = atomIndex_;

    auto inter_residue_bond_count = structure_.bondAtomList.size() / 2;
    while (interBondIndex_ < inter_residue_bond_count) {
        auto atom1 = static_cast<size_t>(structure_.bondAtomList[interBondIndex_ * 2 + 0]);
        auto atom2 = static_cast<size_t>(structure_.bondAtomList[interBondIndex_ * 2 + 1]);

        // We are below the atoms we care about
        if ((atom1 < atomSkip_) || (atom2 < atomSkip_)) {
            interBondIndex_++;
            continue;
        }

        break;
    }

    read(frame);
}

void MMTFFormat::read(Frame& frame) {
    const auto& cell = structure_.unitCell;
    if (structure_.unitCell.size() == 6) {
        Vector3D lengths = {static_cast<double>(cell[0]), static_cast<double>(cell[1]), static_cast<double>(cell[2])};
        Vector3D angles = {static_cast<double>(cell[3]), static_cast<double>(cell[4]), static_cast<double>(cell[5])};
        frame.set_cell({lengths, angles});
    }

    if (!mmtf::isDefaultValue(structure_.title)) {
        frame.set("name", structure_.title);
    }

    if (!mmtf::isDefaultValue(structure_.structureId)) {
        frame.set("pdb_idcode", structure_.structureId);
    }

    if (!mmtf::isDefaultValue(structure_.depositionDate)) {
        frame.set("deposition_date", structure_.depositionDate);
    }

    read_model(frame);
    apply_symmetry(frame);

    atomSkip_ = atomIndex_;
}

void MMTFFormat::read_model(Frame& frame) {
    auto modelChainCount = static_cast<size_t>(structure_.chainsPerModel[modelIndex_]);

    // The following block of code simply counts the total number of atoms in the model
    size_t natoms = 0;
    auto chain_index = chainIndex_;
    auto group_index = groupIndex_;
    for (size_t j = 0; j < modelChainCount; j++) {
        auto chainGroupCount = static_cast<size_t>(structure_.groupsPerChain[chain_index]);
        for (size_t k = 0; k < chainGroupCount; k++) {
            auto groupType = static_cast<size_t>(structure_.groupTypeList[group_index]);
            const auto& group = structure_.groupList[groupType];
            natoms += group.atomNameList.size();
            group_index++;
        }
        chain_index++;
    }

    frame.resize(natoms);
    auto positions = frame.positions();

    // Read the structure iterating over the chains in the model, then the
    // residues/groups in the chain and finally the atoms in the residue/group
    for (size_t j = 0; j < modelChainCount; j++) {

        auto current_assembly = find_assembly();

        auto chainGroupCount = static_cast<size_t>(structure_.groupsPerChain[chainIndex_]);
        // A group is like a residue or other molecule in a PDB file.
        for (size_t k = 0; k < chainGroupCount; k++) {
            // groupType is a 'pointer' to the appropriate residue information in the structure.
            auto group_type = static_cast<size_t>(structure_.groupTypeList[groupIndex_]);
            auto residue = create_residue(current_assembly, group_type);

            read_group(frame, group_type, residue, positions);
            frame.add_residue(std::move(residue));

            add_inter_residue_bonds(frame);

            groupIndex_++;
        }
        chainIndex_++;
    }
    modelIndex_++;
}

std::string MMTFFormat::find_assembly() {
    // Unfortunately we must loop through the assembly lists to find which
    // one our current chain belongs to. Fortunately, these lists are fairly
    // short in the vast majority of cases.

    for (const auto& assembly : structure_.bioAssemblyList) {
        for (const auto& transform : assembly.transformList) {
            for (auto id : transform.chainIndexList) {
                if (static_cast<size_t>(id) == chainIndex_) {
                    return "bio" + assembly.name;
                }
            }
        }
    }

    return "";
}

Residue MMTFFormat::create_residue(const std::string& current_assembly, size_t group_type) {
    const auto& group = structure_.groupList[group_type];

    auto groupId = static_cast<int64_t>(structure_.groupIdList[groupIndex_]);
    auto residue = Residue(group.groupName, groupId);
    residue.set("composition_type", group.chemCompType);
    residue.set("is_standard_pdb", !mmtf::is_hetatm(group.chemCompType.c_str()));

    if (groupIndex_ < structure_.secStructList.size()) {
        set_secondary(residue, structure_.secStructList[groupIndex_]);
    }

    // If the name of the current assembly is defined in the MMTF file.
    // Bioassemblies are optional, however.
    if (!current_assembly.empty()) {
        residue.set("assembly", current_assembly);
    }

    // This is a string in MMTF, differs from the name as then increments linearly
    // For example, the fourth chainid in ( A B A B) would be D, not B (the chainname)
    residue.set("chainid", structure_.chainIdList[chainIndex_]);
    residue.set("chainindex", chainIndex_);
    if (!mmtf::isDefaultValue(structure_.chainNameList)) {
        residue.set("chainname", structure_.chainNameList[chainIndex_]);
    }

    return residue;
}

void MMTFFormat::read_group(Frame& frame, size_t group_type, Residue& residue, span<Vector3D> positions) {

    const auto& group = structure_.groupList[group_type];

    // map between the index of the atom in this group and the global
    // index in the frame. This is needed to add bonds from this group
    std::vector<size_t> global_indexes;
    global_indexes.reserve(group.atomNameList.size());

    for (size_t l = 0; l < group.atomNameList.size(); l++) {
        auto atom = Atom(group.atomNameList[l], group.elementList[l]);
        atom.set_charge(static_cast<double>(group.formalChargeList[l]));

        const auto& altLocList = structure_.altLocList;
        if (!mmtf::isDefaultValue(altLocList)
            && altLocList[atomIndex_] != ' '
            && altLocList[atomIndex_] != 0x00
        ) {
            atom.set("altloc", std::string(1, altLocList[atomIndex_]));
        }

        auto id = atom_id(atomIndex_);
        global_indexes.emplace_back(id);
        frame[id] = std::move(atom);
        residue.add_atom(id);

        positions[id][0] = static_cast<double>(structure_.xCoordList[atomIndex_]);
        positions[id][1] = static_cast<double>(structure_.yCoordList[atomIndex_]);
        positions[id][2] = static_cast<double>(structure_.zCoordList[atomIndex_]);

        atomIndex_++;
    }

    for (size_t l = 0; l < group.bondOrderList.size(); l++) {
        // atom1 and atom2 refer to data in group, so we need to
        // translate them to the global atomic index
        auto atom1 = static_cast<size_t>(group.bondAtomList[l * 2]);
        auto atom2 = static_cast<size_t>(group.bondAtomList[l * 2 + 1]);

        frame.add_bond(
            global_indexes[atom1],
            global_indexes[atom2],
            bond_order_to_chemfiles(group.bondOrderList[l])
        );
    }
}

void MMTFFormat::add_inter_residue_bonds(Frame& frame) {
    auto inter_residue_bond_count = structure_.bondAtomList.size() / 2;

    // Add additional global (not by group) bonds
    while (interBondIndex_ < inter_residue_bond_count) {
        auto atom1 = static_cast<size_t>(structure_.bondAtomList[interBondIndex_ * 2]);
        auto atom2 = static_cast<size_t>(structure_.bondAtomList[interBondIndex_ * 2 + 1]);

        // We are above the atoms we care about
        if (atom1 > atomIndex_ || atom2 > atomIndex_) {
            break;
        }

        frame.add_bond(atom_id(atom1), atom_id(atom2));
        interBondIndex_++;
    }
}

void MMTFFormat::apply_symmetry(Frame& frame) {
    const auto original_size = frame.size();
    const auto original_bond_size = frame.topology().bonds().size();

    using bond_w_order = std::pair<Bond, Bond::BondOrder>;
    std::vector<bond_w_order> bonds_to_add;

    for (const auto& assembly : structure_.bioAssemblyList) {

        for (const auto& transform : assembly.transformList) {

            std::unordered_set<double> chains_to_transform;
            for (auto id : transform.chainIndexList) {
                chains_to_transform.insert(static_cast<double>(id));
            }

            // ncs is a 4x4 matrix stored in column major order.
            const auto& ncs = transform.matrix;
            auto rotation = Matrix3D(
                static_cast<double>(ncs[0]), static_cast<double>(ncs[4]), static_cast<double>(ncs[8]),
                static_cast<double>(ncs[1]), static_cast<double>(ncs[5]), static_cast<double>(ncs[9]),
                static_cast<double>(ncs[2]), static_cast<double>(ncs[6]), static_cast<double>(ncs[10])
            );
            auto translation = Vector3D(
                static_cast<double>(ncs[3]),
                static_cast<double>(ncs[7]),
                static_cast<double>(ncs[11])
            );

            if (rotation == Matrix3D::unit() && translation == Vector3D()) {
                continue;
            }

            rotation = rotation.invert();

            // Makes an atom into its symmetry partner
            std::vector<size_t> old_to_sym(original_size, 0);

            std::vector<Residue> residues_to_add;
            for (const auto& residue : frame.topology().residues()) {
                auto assembly_s = residue.get("assembly");

                if (!assembly_s || assembly_s->as_string() != "bio" + assembly.name) {
                    continue;
                }

                auto chainindex = residue.get("chainindex");
                if (!chainindex || chains_to_transform.count(chainindex->as_double()) == 0) {
                    continue;
                }

                // Copy over everything except the atoms
                auto new_residue = Residue(residue.name(), *residue.id());
                for (const auto& prop : residue.properties()) {
                    new_residue.set(prop.first, prop.second);
                }

                // Avoid using this chain in future symmetry operations
                new_residue.set("chainindex", -1.0);

                for (auto atom_id : residue) {

                    // Ensure that the current atom is not a result of a symmetry operation
                    if (atom_id >= original_size) {
                        continue;
                    }

                    auto new_atom = frame[atom_id];
                    auto new_position = rotation * frame.positions()[atom_id] + translation;

                    frame.add_atom(std::move(new_atom), new_position);
                    new_residue.add_atom(frame.size() - 1);
                    old_to_sym[atom_id] = frame.size() - 1;
                }

                residues_to_add.emplace_back(new_residue);
            }

            for (auto&& residue : residues_to_add) {
                frame.add_residue(std::move(residue));
            }

            for (size_t i = 0; i < original_bond_size; ++i) {
                const auto& bond = frame.topology().bonds()[i];

                // bonds should be sorted so that when we hit original size, we're done
                if (bond[0] >= original_size || bond[1] >= original_size) {
                    break;
                }

                auto new_bond_0 = old_to_sym[bond[0]];
                auto new_bond_1 = old_to_sym[bond[1]];

                // zero simply means that the atom has no partner
                if (new_bond_0 == 0 || new_bond_1 == 0) {
                    continue;
                }

                bonds_to_add.push_back({ {new_bond_0, new_bond_1}, frame.topology().bond_orders()[i] });
            }
        }
    }

    for (auto& bond : bonds_to_add) {
        frame.add_bond(bond.first[0], bond.first[1], bond.second);
    }
}

void MMTFFormat::write(const Frame& frame) {
    structure_.numModels++;
    structure_.chainsPerModel.emplace_back(0);
    structure_.numAtoms += static_cast<int32_t>(frame.size());

    if (mmtf::isDefaultValue(structure_.unitCell)) {
        const auto& cell = frame.cell();
        auto lengths = cell.lengths();
        auto angles = cell.angles();
        structure_.unitCell.resize(6);
        structure_.unitCell[0] = static_cast<float>(lengths[0]);
        structure_.unitCell[1] = static_cast<float>(lengths[1]);
        structure_.unitCell[2] = static_cast<float>(lengths[2]);
        structure_.unitCell[3] = static_cast<float>(angles[0]);
        structure_.unitCell[4] = static_cast<float>(angles[1]);
        structure_.unitCell[5] = static_cast<float>(angles[2]);

        unitcellForWrite_ = cell;
    } else if (unitcellForWrite_ != frame.cell()) {
        warning("MMTF Writer", "the MMTF format only allows one unit cell to be defined for all models, using the first one defined");
    }

    const auto& topology = frame.topology();

    // pre-allocate some memory
    auto size = static_cast<size_t>(structure_.numAtoms);
    structure_.xCoordList.reserve(size);
    structure_.yCoordList.reserve(size);
    structure_.zCoordList.reserve(size);
    structure_.atomIdList.reserve(size);

    structure_.groupList.reserve(structure_.groupList.size() + topology.residues().size());

    new_atom_indexes_.clear();
    new_atom_indexes_.resize(frame.size(), -1);

    std::string previous_chainId;
    std::string previous_chainName;
    bool got_at_least_one_chain_in_this_frame = false;
    for (const auto& residue: topology.residues()) {
        // WARNING: this assumes that residues in the same chain are contiguous
        //
        // [residue[chainid="A"], residue[chainid="B"], residue[chainid="A"]]
        // will result in three chains
        const auto& chainId = residue.get<Property::STRING>("chainid").value_or("");
        const auto& chainName = residue.get<Property::STRING>("chainname").value_or("");
        if (chainName != previous_chainName ||
            chainId != previous_chainId ||
            !got_at_least_one_chain_in_this_frame
        ) {
            got_at_least_one_chain_in_this_frame = true;
            structure_.numChains++;
            structure_.chainsPerModel.back()++;

            // No residues in this chain yet
            structure_.groupsPerChain.emplace_back(0);

            structure_.chainIdList.emplace_back(chainId);
            structure_.chainNameList.emplace_back(chainName);

            previous_chainId = chainId;
            previous_chainName = chainName;
        }

        this->add_residue_to_structure(frame, residue);
    }

    // For all atoms not inside a residue, add them inside a new residue,
    // containing only this atom
    auto not_in_residue = std::vector<size_t>();
    for (size_t i=0; i<frame.size(); i++) {
        if (!topology.residue_for_atom(i)) {
            not_in_residue.emplace_back(i);
        }
    }

    if (!not_in_residue.empty()) {
        // add all of these new articficial residues in a separate chain
        structure_.numChains++;
        structure_.chainsPerModel.back()++;
        structure_.groupsPerChain.emplace_back(0);
        structure_.chainIdList.emplace_back("");
        structure_.chainNameList.emplace_back("");

        for (auto i: not_in_residue) {
            auto residue = Residue("");
            residue.add_atom(i);
            this->add_residue_to_structure(frame, residue);
        }
    }

    mmtf::BondAdder add_mmtf_bond(structure_);
    const auto& bonds = topology.bonds();
    const auto& bond_orders = topology.bond_orders();
    for (size_t i = 0; i < bonds.size(); ++i) {
        add_mmtf_bond(
            new_atom_indexes_[bonds[i][0]],
            new_atom_indexes_[bonds[i][1]],
            bond_order_to_mmtf(bond_orders[i])
        );
    }

    atomSkip_ += frame.size();
}

void MMTFFormat::add_residue_to_structure(const Frame& frame, const Residue& residue) {

    structure_.numGroups++;
    structure_.groupsPerChain.back() += 1;

    auto groupType = static_cast<int32_t>(structure_.groupList.size());
    structure_.groupTypeList.emplace_back(groupType);

    int32_t groupId = residue.id() ? static_cast<int32_t>(residue.id().value()) : -1;
    structure_.groupIdList.emplace_back(groupId);

    auto group = mmtf::GroupType();
    group.groupName = residue.name();

    auto composition_type = residue.get<Property::STRING>("composition_type").value_or("other");
    group.chemCompType = std::move(composition_type);

    // pre-allocate some memory
    group.formalChargeList.reserve(residue.size());
    group.atomNameList.reserve(residue.size());
    group.elementList.reserve(residue.size());

    const auto& positions = frame.positions();
    for (auto i: residue) {
        const auto& atom = frame[i];
        group.formalChargeList.emplace_back(atom.charge());

        if (atom.name().size() > 5) {
            warning("MMTF Writer",
                "atom name '{}' is too long for MMTF format, it will be truncated",
                atom.name()
            );
        }
        group.atomNameList.emplace_back(atom.name().substr(0, 5));

        if (atom.type().size() > 3) {
            warning("MMTF Writer",
                "atom type '{}' is too long for MMTF format, it will be truncated",
                atom.type()
            );
        }
        group.elementList.emplace_back(atom.type().substr(0, 3));

        new_atom_indexes_[i] = static_cast<int32_t>(structure_.xCoordList.size());
        structure_.atomIdList.emplace_back(atomSkip_ + i + 1);
        structure_.xCoordList.emplace_back(positions[i][0]);
        structure_.yCoordList.emplace_back(positions[i][1]);
        structure_.zCoordList.emplace_back(positions[i][2]);
    }
    structure_.groupList.emplace_back(std::move(group));
}

// A function to translate from the index in MMTF lists to an atom id
// suitable for chemfiles: starts at 0 for each model, and correspond to
// the initial atom index if it exists.
size_t MMTFFormat::atom_id(size_t mmtf_id) {
    if (!mmtf::isDefaultValue(structure_.atomIdList)) {
        assert(mmtf_id < structure_.atomIdList.size());
        assert(structure_.atomIdList[mmtf_id] > 0);
        size_t id;
        if (new_atom_indexes_.empty()) {
            // atom ids are well-behaved, no reordering necessary
            id = static_cast<size_t>(structure_.atomIdList[mmtf_id]) - 1;
        } else {
            assert(structure_.atomIdList.size() == new_atom_indexes_.size());
            auto listedId = structure_.atomIdList[mmtf_id];
            // use the fact that indexes are sorted for faster search
            auto it =
                std::lower_bound(new_atom_indexes_.begin(), new_atom_indexes_.end(), listedId);
            // it is guaranteed that the lower bound is an exact match
            assert(it != new_atom_indexes_.end() && *it == listedId);
            id = static_cast<size_t>(std::distance(new_atom_indexes_.begin(), it));
        }
        assert(atomSkip_ <= id);
        return id - atomSkip_;
    } else {
        assert(atomSkip_ <= mmtf_id);
        return mmtf_id - atomSkip_;
    }
}

int8_t bond_order_to_mmtf(Bond::BondOrder order) {
    auto warn_bad_bond_order = [](const char* name) {
        warning("MMTF Writer",
            "bond order '{}' can not be represented in MMTF, defaulting to single bond",
            name
        );
        return static_cast<int8_t>(1);
    };

    switch(order) {
    case Bond::BondOrder::SINGLE:
        return 1;
    case Bond::BondOrder::DOUBLE:
        return 2;
    case Bond::BondOrder::TRIPLE:
        return 3;
    case Bond::BondOrder::QUADRUPLE:
        return 4;
    case Bond::BondOrder::UNKNOWN:
        return -1;
    case Bond::BondOrder::QUINTUPLET:
        return warn_bad_bond_order("quintuplet");
    case Bond::BondOrder::AMIDE:
        return warn_bad_bond_order("amide");
    case Bond::BondOrder::AROMATIC:
        return warn_bad_bond_order("aromatic");
    case Bond::BondOrder::UP:
        return warn_bad_bond_order("up");
    case Bond::BondOrder::DOWN:
        return warn_bad_bond_order("down");
    case Bond::BondOrder::DATIVE_L:
        return warn_bad_bond_order("dative_l");
    case Bond::BondOrder::DATIVE_R:
        return warn_bad_bond_order("dative_r");
    default:
        unreachable();
    }
}

Bond::BondOrder bond_order_to_chemfiles(int32_t order) {
    switch(order) {
    case 1:
        return Bond::SINGLE;
    case 2:
        return Bond::DOUBLE;
    case 3:
        return Bond::TRIPLE;
    case 4:
        return Bond::QUADRUPLE;
    case -1:
        return Bond::UNKNOWN;
    default:
        warning("MMTF Reader", "unexpected bond order from MMTF '{}'", order);
        return Bond::UNKNOWN;
    }
}

void set_secondary(Residue& residue, int32_t code) {
    // We use the raw values here to directly match the MMTF spec
    // See https://github.com/rcsb/mmtf/blob/master/spec.md#secstructlist
    switch(code) {
    case 0:
        residue.set("secondary_structure", "pi helix");
        break;
    case 1:
        residue.set("secondary_structure", "bend");
        break;
    case 2:
        residue.set("secondary_structure", "alpha helix");
        break;
    case 3:
        residue.set("secondary_structure", "extended");
        break;
    case 4:
        residue.set("secondary_structure", "3-10 helix");
        break;
    case 5:
        residue.set("secondary_structure", "bridge");
        break;
    case 6:
        residue.set("secondary_structure", "turn");
        break;
    case 7:
        residue.set("secondary_structure", "coil");
        break;
    case -1:
        break;
    default:
        warning("MMTF Reader", "unknown secondary structure code '{}'", code);
        break;
    }
}
