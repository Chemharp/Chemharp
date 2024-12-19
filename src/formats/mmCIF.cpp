// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cassert>
#include <cstddef>
#include <cstdint>

#include <map>
#include <array>
#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include <string_view>

#include "chemfiles/types.hpp"
#include "chemfiles/utils.hpp"
#include "chemfiles/parse.hpp"
#include "chemfiles/error_fmt.hpp"
#include "chemfiles/external/optional.hpp"

#include "chemfiles/File.hpp"
#include "chemfiles/Atom.hpp"
#include "chemfiles/Frame.hpp"
#include "chemfiles/Property.hpp"
#include "chemfiles/Residue.hpp"
#include "chemfiles/Topology.hpp"
#include "chemfiles/UnitCell.hpp"
#include "chemfiles/FormatMetadata.hpp"

// WARNING UGLY HACK!
#include "chemfiles/formats/PDB.hpp"

#include "chemfiles/formats/mmCIF.hpp"

using namespace chemfiles;

template<> const FormatMetadata& chemfiles::format_metadata<mmCIFFormat>() {
    static FormatMetadata metadata;
    metadata.name = "mmCIF";
    metadata.extension = ".mmcif";
    metadata.description = "Crystallographic Information Framework files for MacroMolecules";
    metadata.reference = "http://mmcif.wwpdb.org/";

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

/// CIF files store which digits are insignificant inside parenthesis, this
/// function handle this
static double cif_to_double(std::string_view line);

void mmCIFFormat::init_() {
    if (file_.mode() == File::WRITE) {
        return;
    }

    if (file_.mode() == File::APPEND) {
        throw file_error("cannot open mmCIF files in append ('a') mode");
    }

    Vector3D lengths;
    Vector3D angles = {90, 90, 90};

    bool in_loop = false;
    size_t current_index = 0;
    while (!file_.eof()) {
        auto line = file_.readline();

        if (line.find("loop_") != std::string::npos) {
            in_loop = true;
            continue;
        }

        if (line.empty() || line[0] == '#') {
            in_loop = false;
            continue;
        }

        auto line_split = split(line, ' ');
        if (line_split.size() > 1 && line[0] == '_') {
            in_loop = false;
        }

        if (line_split[0] == "_cell_length_a" || line_split[0] == "_cell.length_a") {
            lengths[0] = cif_to_double(line_split[1]);
        }

        if (line_split[0] == "_cell_length_b" || line_split[0] == "_cell.length_b") {
            lengths[1] = cif_to_double(line_split[1]);
        }

        if (line_split[0] == "_cell_length_c" || line_split[0] == "_cell.length_c") {
            lengths[2] = cif_to_double(line_split[1]);
        }

        if (line_split[0] == "_cell_angle_alpha" || line_split[0] == "_cell.angle_alpha") {
            angles[0] = cif_to_double(line_split[1]);
        }

        if (line_split[0] == "_cell_angle_beta" || line_split[0] == "_cell.angle_beta") {
            angles[1] = cif_to_double(line_split[1]);
        }

        if (line_split[0] == "_cell_angle_gamma" || line_split[0] == "_cell.angle_gamma") {
            angles[2] = cif_to_double(line_split[1]);
        }

        if (line_split[0] == "_entry.id") {
            pdb_idcode_ = std::string(line_split[1]);
        }

        if (line_split[0] == "_struct.title") {
            auto trimmed = trim(line.substr(13));
            name_ = trimmed.size() > 2 ? std::string(trimmed.substr(1, trimmed.size() - 2)) : "";
        }

        if (in_loop && line_split[0].find("_atom_site.") != std::string::npos) {
            auto atom_label = std::string(line_split[0].substr(11));
            atom_site_map_[atom_label] = current_index++;
            break;
        }
    }

    if (current_index == 0) {
        throw format_error("could not find _atom_site loop in '{}'", file_.path());
    }

    cell_ = UnitCell(lengths, angles);

    auto position = file_.tellpos();
    auto line = file_.readline();

    do {
        if (line.find("_atom_site") != std::string::npos) {
            auto atom_label = std::string(trim(line).substr(11));
            atom_site_map_[atom_label] = current_index++;

            position = file_.tellpos();
            line = file_.readline();

            continue;
        }

        // Anything else ends the block
        break;
    } while (!file_.eof());

    // After this block ends, we have the start of coordinates
    steps_positions_.push_back(position);

    if (atom_site_map_.find("type_symbol") == atom_site_map_.end()) {
        throw format_error("could not find _atom_site.type_symbol in '{}'", file_.path());
    }

    if (atom_site_map_.find("Cartn_x") == atom_site_map_.end()) {
        throw format_error("could not find _atom_site.Cartn_x in '{}'", file_.path());
    }

    auto model_position = atom_site_map_.find("pdbx_PDB_model_num");
    // Do we have a special extension for multiple modes?
    if (model_position == atom_site_map_.end()) {
        // If not, we are done
        file_.seekpos(steps_positions_[0]);
        return;
    }

    // Ok, let's look at the sites now to note where models start
    auto last_position = parse<size_t>(split(line, ' ')[model_position->second]);

    do {
        position = file_.tellpos();
        line = file_.readline();

        // a break in the text ends the models
        if (line.empty() || line == "line_" || line[0] == '#') {
            break;
        }

        auto line_split = split(line, ' ');
        size_t current_position = parse<size_t>(line_split[model_position->second]);

        if (current_position != last_position) {
            steps_positions_.push_back(position);
            last_position = current_position;
        }
    } while (!file_.eof());

    file_.seekpos(steps_positions_[0]);
}

size_t mmCIFFormat::nsteps() {
    return steps_positions_.size();
}

void mmCIFFormat::read_step(const size_t step, Frame& frame) {
    assert(step < steps_positions_.size());
    file_.seekpos(steps_positions_[step]);
    read(frame);
}

void mmCIFFormat::read(Frame& frame) {
    map_residues_indexes.clear();
    residues_.clear();
    frame.set_cell(cell_);

    if (!name_.empty()) {
        frame.set("name", name_);
    }

    if (!pdb_idcode_.empty()) {
        frame.set("pdb_idcode", pdb_idcode_);
    }

    // The following map operations can be moved to the constructor
    // and stored with optional. I don't know which solution is better.

    // These are required for atoms
    auto type_symbol = atom_site_map_.at("type_symbol");

    // This has two names...
    auto label_atom_id = atom_site_map_.find("label_atom_id");
    if (label_atom_id == atom_site_map_.end()) {
        label_atom_id = atom_site_map_.find("label");
    }

    // Other atom properties
    auto group_pdb = atom_site_map_.find("group_PDB");
    auto label_alt_id = atom_site_map_.find("label_alt_id");
    auto formal_charge = atom_site_map_.find("formal_charge");

    // Positions
    auto cartn_x = atom_site_map_.at("Cartn_x");
    auto cartn_y = atom_site_map_.at("Cartn_y");
    auto cartn_z = atom_site_map_.at("Cartn_z");

    // Residue properties
    auto label_comp_id = atom_site_map_.find("label_comp_id");
    auto label_asym_id = atom_site_map_.find("label_asym_id");
    auto auth_asym_id = atom_site_map_.find("auth_asym_id");
    auto label_seq_id = atom_site_map_.find("label_seq_id");
    auto label_entity_id = atom_site_map_.find("label_entity_id");

    auto model_position = atom_site_map_.find("pdbx_PDB_model_num");

    auto position = file_.tellpos();

    size_t last_position = 0;
    if (model_position != atom_site_map_.end()) {
        auto line = file_.readline();
        last_position = parse<size_t>(split(line, ' ')[model_position->second]);
        // Reset file position so that the loop below can start by reading the
        // first line
        file_.seekpos(position);
    }

    while (!file_.eof()) {
        auto line = file_.readline();
        auto line_split = split(line, ' ');
        if (line.empty() || line == "loop_" || line[0] == '#') {
            break;
        }

        if (line_split.size() != atom_site_map_.size()) {
            throw format_error("line '{}' has {} items not {}",
                line, line_split.size(), atom_site_map_.size()
            );
        }

        size_t current_position = 0;
        if (model_position != atom_site_map_.end()) {
            current_position = parse<size_t>(line_split[model_position->second]);
        }

        if (current_position != last_position) {
            break;
        }

        auto atom = Atom(
            std::string(line_split[label_atom_id->second]),
            std::string(line_split[type_symbol])
        );

        if (label_alt_id != atom_site_map_.end() &&
            line_split[label_alt_id->second] != ".") {
            atom.set("altloc", std::string(line_split[label_alt_id->second]));
        }

        if (formal_charge != atom_site_map_.end()) {
            atom.set_charge(cif_to_double(line_split[formal_charge->second]));
        }

        auto x = cif_to_double(line_split[cartn_x]);
        auto y = cif_to_double(line_split[cartn_y]);
        auto z = cif_to_double(line_split[cartn_z]);
        frame.add_atom(std::move(atom), Vector3D(x, y, z));

        position = file_.tellpos();

        if (label_comp_id == atom_site_map_.end() || label_asym_id == atom_site_map_.end()) {
            continue;
        }

        auto atom_id = frame.size() - 1;
        int64_t resid = 0;
        auto resid_text = line_split[label_seq_id->second];

        try {
            if (resid_text == ".") { // In this case, we need to use the entity id
                resid = parse<int64_t>(line_split[label_entity_id->second]);
            } else {
                resid = parse<int64_t>(line_split[label_seq_id->second]);
            }
        } catch (const Error& e) {
            throw format_error("invalid CIF residue or entity numeric: {}", e.what());
        }

        auto chainid = std::string(line_split[label_asym_id->second]);

        if (map_residues_indexes.find({chainid, resid}) == map_residues_indexes.end()) {

            auto name = line_split[label_comp_id->second];
            Residue residue(std::string(name), resid);
            residue.add_atom(atom_id);

            // This will be saved as a string on purpose to match MMTF
            if (label_asym_id != atom_site_map_.end()) {
                residue.set("chainid", chainid);
            }

            if (auth_asym_id != atom_site_map_.end()) {
                residue.set("chainname", std::string(line_split[auth_asym_id->second]));
            }

            if (group_pdb != atom_site_map_.end()) {
                residue.set("is_standard_pdb", line_split[group_pdb->second] == "ATOM");
            }

            map_residues_indexes.emplace(std::make_pair(chainid, resid), residues_.size());
            residues_.emplace_back(std::move(residue));
        } else {
            // Just add this atom to the residue
            residues_[map_residues_indexes.at({ chainid, resid })].add_atom(atom_id);
        }
    }

    // Reset state to previous line
    file_.seekpos(position);

    for (const auto& residue: residues_) {
        frame.add_residue(residue);
    }

    // Only link if we are reading mmCIF
    if (model_position != atom_site_map_.end()) {
        // Cross format talk! Forgive me!
        PDBFormat::link_standard_residue_bonds(frame);
    }
}

void mmCIFFormat::write(const Frame& frame) {
    if (models_ == 0) {
        file_.print("# generated by Chemfiles\n");
        file_.print("#\n");
        auto lengths = frame.cell().lengths();
        file_.print("_cell.length_a {:#g}\n", lengths[0]);
        file_.print("_cell.length_b {:#g}\n", lengths[1]);
        file_.print("_cell.length_c {:#g}\n", lengths[2]);
        auto angles = frame.cell().angles();
        file_.print("_cell.length_alpha {:#g}\n", angles[0]);
        file_.print("_cell.length_beta  {:#g}\n", angles[1]);
        file_.print("_cell.length_gamma {:#g}\n", angles[2]);
        file_.print("#\n");
        file_.print("loop_\n");
        file_.print("_atom_site.group_PDB\n");
        file_.print("_atom_site.id\n");
        file_.print("_atom_site.type_symbol\n");
        file_.print("_atom_site.label_atom_id\n");
        file_.print("_atom_site.label_alt_id\n");
        file_.print("_atom_site.label_comp_id\n");
        file_.print("_atom_site.label_asym_id\n");
        file_.print("_atom_site.label_seq_id\n");
        file_.print("_atom_site.Cartn_x\n");
        file_.print("_atom_site.Cartn_y\n");
        file_.print("_atom_site.Cartn_z\n");
        file_.print("_atom_site.pdbx_formal_charge\n");
        file_.print("_atom_site.auth_asym_id\n");
        file_.print("_atom_site.pdbx_PDB_model_num\n");
    }

    models_++;

    const auto& topology = frame.topology();
    const auto& positions = frame.positions();
    for (size_t i = 0; i < frame.size(); ++i) {
        ++atoms_;

        std::string compid = ".";
        std::string asymid = ".";
        std::string seq_id = ".";
        std::string auth_asymid = ".";
        std::string pdbgroup = "HETATM";

        const auto& residue = topology.residue_for_atom(i);
        if (residue) {
            compid = residue->name();

            if (residue->id()) {
                seq_id = std::to_string(*residue->id());
            } else {
                seq_id = "?";
            }

            asymid = residue->get<Property::STRING>("chainid").value_or("?");
            auth_asymid = residue->get<Property::STRING>("chainname").value_or(".");
            if (residue->get<Property::BOOL>("is_standard_pdb").value_or(false)) {
                pdbgroup = "ATOM  ";
            }
        }

        const auto& atom = frame[i];

        file_.print("{} {: <5} {: <2} {: <4} {} {: >3} {} {: >4} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {} {}\n",
                pdbgroup, atoms_, atom.type(), atom.name(), ".", compid,
                asymid, seq_id, positions[i][0], positions[i][1], positions[i][2],
                atom.charge(), auth_asymid, models_
        );
    }

}

double cif_to_double(std::string_view line) {
    // only allocate a new string if we actually find `(` somewhere, otherwise
    // just parse it
    std::string line_string;

    const auto* open = std::find(line.begin(), line.end(), '(');
    if (open != line.end()) {
        line_string = std::string(line);
        line_string.erase(std::remove(line_string.begin(), line_string.end(), '('), line_string.end());
        line_string.erase(std::remove(line_string.begin(), line_string.end(), ')'), line_string.end());

        line = line_string;
    }

    return parse<double>(line);
}
