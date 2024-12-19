// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cmath>
#include <cassert>
#include <cstddef>
#include <cstdint>

#include <map>
#include <array>
#include <deque>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <string_view>

#include "chemfiles/types.hpp"
#include "chemfiles/utils.hpp"
#include "chemfiles/parse.hpp"
#include "chemfiles/warnings.hpp"
#include "chemfiles/error_fmt.hpp"
#include "chemfiles/external/optional.hpp"

#include "chemfiles/File.hpp"
#include "chemfiles/Atom.hpp"
#include "chemfiles/Frame.hpp"
#include "chemfiles/Property.hpp"
#include "chemfiles/Residue.hpp"
#include "chemfiles/Topology.hpp"
#include "chemfiles/UnitCell.hpp"
#include "chemfiles/Connectivity.hpp"
#include "chemfiles/FormatMetadata.hpp"

#include "chemfiles/formats/PDB.hpp"
#include "chemfiles/pdb_connectivity.hpp"

using namespace chemfiles;

template<> const FormatMetadata& chemfiles::format_metadata<PDBFormat>() {
    static FormatMetadata metadata;
    metadata.name = "PDB";
    metadata.extension = ".pdb";
    metadata.description = "PDB (RCSB Protein Data Bank) text format";
    metadata.reference = "http://www.wwpdb.org/documentation/file-format";

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

bool chemfiles::operator==(const FullResidueId& lhs, const FullResidueId& rhs) {
    return lhs.chain == rhs.chain &&
           lhs.resid == rhs.resid &&
           lhs.resname == rhs.resname &&
           lhs.insertion_code == rhs.insertion_code;
}

bool chemfiles::operator<(const FullResidueId& lhs, const FullResidueId& rhs) {
    if (lhs.chain != rhs.chain) {
        return lhs.chain < rhs.chain;
    } else if (lhs.resid != rhs.resid) {
        return lhs.resid < rhs.resid;
    } else if (lhs.insertion_code != rhs.insertion_code) {
        return lhs.insertion_code < rhs.insertion_code;
    } else if (lhs.resname != rhs.resname) {
        // Check resname last, since resid should already have taken care of
        // identical residues. Files containing residues with different names
        // but the same resid can be created when multiple PDB are merged
        // together
        return lhs.resname < rhs.resname;
    } else {
        return false;
    }
}

/// Check the number of digits before the decimal separator to be sure than we
/// can represent them. In case of error, use the given `context` in the error
/// message
static void check_values_size(const Vector3D& values, unsigned width, const std::string& context);

// PDB record handled by chemfiles. Any record not in this enum are not yet
// implemented.
enum class Record {
    // Records containing summary data
    HEADER,
    TITLE,
    // Records containing usefull data
    CRYST1,
    ATOM,
    HETATM,
    CONECT,
    // Beggining of model.
    MODEL,
    // End of model.
    ENDMDL,
    // End of chain. May increase atom count
    TER,
    // End of file
    END,
    // Secondary structure
    HELIX,
    SHEET,
    TURN,
    // Ignored records
    IGNORED_,
    // Unknown record type
    UNKNOWN_,
};

// Get the record type for a line.
static Record get_record(std::string_view line);

void PDBFormat::read_next(Frame& frame) {
    residues_.clear();
    atom_offsets_.clear();

    uint64_t position;
    bool got_end = false;
    while (!got_end && !file_.eof()) {
        auto line = file_.readline();
        auto record = get_record(line);
        std::string name;
        switch (record) {
        case Record::HEADER:
            if (line.size() >= 50) {
                frame.set("classification", std::string(trim(line.substr(10, 40))));
            }
            if (line.size() >= 59) {
                frame.set("deposition_date", std::string(trim(line.substr(50, 9))));
            }
            if (line.size() >= 66) {
                frame.set("pdb_idcode", std::string(trim(line.substr(62, 4))));
            }
            continue;
        case Record::TITLE:
            if (line.size() < 11) {continue;}
            // get previous frame name (from a previous TITLE record) and
            // append to it
            name = frame.get<Property::STRING>("name").value_or("");
            name = name.empty() ? "" : name + " ";
            frame.set("name", name + std::string(trim(line.substr(10 , 70))));
            continue;
        case Record::CRYST1:
            read_CRYST1(frame, line);
            continue;
        case Record::ATOM:
            read_ATOM(frame, line, false);
            continue;
        case Record::HETATM:
            read_ATOM(frame, line, true);
            continue;
        case Record::CONECT:
            read_CONECT(frame, line);
            continue;
        case Record::MODEL:
            models_++;
            continue;
        case Record::ENDMDL:
            // Check if the next record is an `END` record
            if (!file_.eof()) {
                position = file_.tellpos();
                line = file_.readline();
                file_.seekpos(position);
                if (get_record(line) == Record::END) {
                    // If this is the case, wait for this next record
                    continue;
                }
            }
            // Else we have read a frame
            got_end = true;
            continue;
        case Record::HELIX:
            read_HELIX(line);
            continue;
        case Record::SHEET:
            read_secondary(line, 17, 28, "SHEET");
            continue;
        case Record::TURN:
            read_secondary(line, 15, 26, "TURN");
            continue;
        case Record::TER:
            if (line.size() >= 12) {
                try {
                    auto ter_serial = decode_hybrid36(5, line.substr(6, 5));
                    if (ter_serial != 0) { // this happens if the ter serial number is blank
                        atom_offsets_.push_back(static_cast<size_t>(ter_serial));
                    }
                } catch (const Error&) {
                    warning("PDB reader", "TER record not numeric: {}", line);
                }
            }
            chain_ended(frame);
            continue;
        case Record::END:
            // We have read a frame!
            got_end = true;
            continue;
        case Record::IGNORED_:
            continue; // Nothing to do
        case Record::UNKNOWN_:
        default:
            if (!file_.eof()) {
                warning("PDB reader", "ignoring unknown record: {}", line);
            }
            continue;
        }
    }

    if (!got_end) {
        warning("PDB reader", "missing END record in file");
    }

    chain_ended(frame);
    link_standard_residue_bonds(frame);
}

void PDBFormat::read_CRYST1(Frame& frame, std::string_view line) {
    assert(line.substr(0, 6) == "CRYST1");
    if (line.length() < 54) {
        throw format_error("CRYST1 record '{}' is too small", line);
    }
    try {
        auto lengths = Vector3D(
            parse<double>(line.substr(6, 9)),
            parse<double>(line.substr(15, 9)),
            parse<double>(line.substr(24, 9))
        );
        auto angles = Vector3D(
            parse<double>(line.substr(33, 7)),
            parse<double>(line.substr(40, 7)),
            parse<double>(line.substr(47, 7))
        );

        frame.set_cell({lengths, angles});
    } catch (const Error&) {
        throw format_error("could not read CRYST1 record '{}'", line);
    }

    if (line.length() >= 55) {
        auto space_group = trim(line.substr(55, 10));
        if (space_group != "P 1" && space_group != "P1") {
            warning("PDB reader", "ignoring custom space group ({}), using P1 instead", space_group);
        }
    }
}

// See http://www.wwpdb.org/documentation/file-format-content/format23/sect5.html
// for definitions of helix types
static const char* HELIX_TYPES[10] = {
    "right-handed alpha helix",
    "right-handed omega helix",
    "right-handed pi helix",
    "right-handed gamma helix",
    "right-handed 3-10 helix",
    "left-handed alpha helix",
    "left-handed omega helix",
    "left-handed gamma helix",
    "2-7 ribbon/helix",
    "polyproline",
};

void PDBFormat::read_HELIX(std::string_view line) {
    if (line.length() < 33 + 5) {
        warning("PDB reader", "HELIX record too short: '{}'", line);
        return;
    }

    auto chain_start = line[19];
    auto chain_end = line[31];
    auto inscode_start = line[25];
    auto inscode_end = line[37];
    auto resname_start = std::string(trim(line.substr(15, 3)));
    auto resname_end = std::string(trim(line.substr(27, 3)));

    int64_t resid_start = 0;
    int64_t resid_end = 0;
    try {
        resid_start = decode_hybrid36(4, line.substr(21, 4));
        resid_end = decode_hybrid36(4, line.substr(33, 4));
    } catch (const Error&) {
        warning("PDB reader", "HELIX record contains invalid numbers: '{}'", line);
        return;
    }

    if (chain_start != chain_end) {
        warning("PDB reader", "HELIX chain {} and {} are not the same", chain_start, chain_end);
        return;
    }

    auto start = FullResidueId {chain_start, resid_start, resname_start, inscode_start};
    auto end = FullResidueId {chain_end, resid_end, resname_end, inscode_end};

    size_t helix_type = 0;
    try {
        helix_type = parse<size_t>(line.substr(38,2));
    } catch (const Error&) {
        warning("PDB reader", "could not parse helix type");
        return;
    }

    if (helix_type <= 10) {
        secinfo_.emplace(start, std::make_pair(end, HELIX_TYPES[helix_type - 1]));
    } else {
        warning("PDB reader", "invalid HELIX type {}", helix_type);
    }
}

void PDBFormat::read_secondary(std::string_view line, size_t i_start, size_t i_end, std::string_view record) {
    if (line.length() < i_end + 10) {
        warning("PDB reader", "secondary structure record too short: '{}'", line);
        return;
    }

    auto resname_start = std::string(trim(line.substr(i_start, 3)));
    auto resname_end = std::string(trim(line.substr(i_end, 3)));

    auto chain_start = line[i_start + 4];
    auto chain_end = line[i_end + 4];

    if (chain_start != chain_end) {
        warning("PDB reader", "{} chain {} and {} are not the same", record, chain_start, chain_end);
        return;
    }

    int64_t resid_start = 0;
    int64_t resid_end = 0;
    try {
        resid_start = decode_hybrid36(4, line.substr(i_start + 5, 4));
        resid_end = decode_hybrid36(4, line.substr(i_end + 5, 4));
    } catch (const Error&) {
        warning("PDB reader",
            "error parsing line: '{}', check {} and {}",
            line, line.substr(i_start + 5, 4), line.substr(i_end + 5, 4)
        );
        return;
    }

    auto inscode_start = line[i_start + 9];
    auto inscode_end = line[i_end + 9];

    auto start = FullResidueId {chain_start, resid_start, resname_start, inscode_start};
    auto end = FullResidueId {chain_end, resid_end, resname_end, inscode_end};

    secinfo_.insert({ start, std::make_pair(end, "extended") });
}

void PDBFormat::read_ATOM(Frame& frame, std::string_view line, bool is_hetatm) {
    assert(line.substr(0, 6) == "ATOM  " || line.substr(0, 6) == "HETATM");

    if (line.length() < 54) {
        throw format_error(
            "{} record is too small: '{}'", line.substr(0, 6), line
        );
    }

    if (atom_offsets_.empty()) {
        try {
            auto initial_offset = decode_hybrid36(5, line.substr(6, 5));

            if (initial_offset <= 0) {
                warning("PDB reader", "{} is too small, assuming id is '1'", initial_offset);
                atom_offsets_.push_back(0);
            } else {
                atom_offsets_.push_back(static_cast<size_t>(initial_offset) - 1);
            }
        } catch(const Error&) {
            warning("PDB reader", "{} is not a valid atom id, assuming '1'", line.substr(6, 5));
			atom_offsets_.push_back(0);
        }
    }

    Atom atom;
    auto name = line.substr(12, 4);
    if (line.length() >= 78) {
        auto type = line.substr(76, 2);
        // Read both atom name and atom type
        atom = Atom(std::string(trim(name)), std::string(trim(type)));
    } else {
        // Read just the atom name and hope for the best.
        atom = Atom(std::string(trim(name)));
    }

    auto altloc = line.substr(16, 1);
    if (altloc != " ") {
        atom.set("altloc", std::string(altloc));
    }

    try {
        auto x = parse<double>(line.substr(30, 8));
        auto y = parse<double>(line.substr(38, 8));
        auto z = parse<double>(line.substr(46, 8));

        frame.add_atom(std::move(atom), Vector3D(x, y, z));
    } catch (const Error&) {
        throw format_error("could not read positions in '{}'", line);
    }

    auto atom_id = frame.size() - 1;
    auto insertion_code = line[26];
    int64_t resid;
    try {
        resid = decode_hybrid36(4, line.substr(22, 4));
    } catch (const Error&) {
        // No residue information
        return;
    }

    auto chain = line[21];
    auto resname = std::string(trim(line.substr(17, 3)));
    auto full_residue_id = FullResidueId {chain, resid, resname, insertion_code};
    if (residues_.count(full_residue_id) == 0) {
        Residue residue(std::move(resname), resid);
        residue.add_atom(atom_id);

        if (insertion_code != ' ') {
            residue.set("insertion_code", std::string(line.substr(26, 1)));
        }

        // Set whether or not the residue is standardized
        residue.set("is_standard_pdb", !is_hetatm);
        // This is saved as a string (instead of a number) on purpose
        // to match MMTF format
        residue.set("chainid", std::string(line.substr(21, 1)));
        // PDB format makes no distinction between chainid and chainname
        residue.set("chainname", std::string(line.substr(21, 1)));

        // segment name is not part of the standard, but something added by
        // CHARM/NAMD in the un-used character range 73-76
        if (line.length() > 72) {
            auto segname = trim(line.substr(72, 4));
            if (segname != "") {
                residue.set("segname", std::string(segname));
            }
        }

        // Are we withing a secondary information sequence?
        if (current_secinfo_) {
            residue.set("secondary_structure", current_secinfo_->second);

            // Are we the end of a secondary information sequence?
            if (current_secinfo_->first == full_residue_id) {
                current_secinfo_ = nullopt;
            }
        }

        // Are we the start of a secondary information sequence?
        auto secinfo_for_residue = secinfo_.find(full_residue_id);
        if (secinfo_for_residue != secinfo_.end()) {
            current_secinfo_ = secinfo_for_residue->second;
            residue.set("secondary_structure", secinfo_for_residue->second.second);
        }

        residues_.emplace(full_residue_id, residue);
    } else {
        // Just add this atom to the residue
        residues_.at(full_residue_id).add_atom(atom_id);
    }
}

void PDBFormat::read_CONECT(Frame& frame, std::string_view line) {
    assert(line.substr(0, 6) == "CONECT");
    auto line_length = trim(line).length();

    // Helper lambdas
    auto add_bond = [&frame, &line](size_t i, size_t j) {
        if (i >= frame.size() || j >= frame.size()) {
            warning("PDB reader",
                "ignoring CONECT ('{}') with atomic indexes bigger than frame size ({})",
                trim(line), frame.size()
            );
            return;
        }
        frame.add_bond(i, j);
    };

    auto read_index = [&line,this](size_t initial) -> size_t {
        try {
            auto pdb_atom_id = decode_hybrid36(5, line.substr(initial, 5));
            auto lower = std::lower_bound(atom_offsets_.begin(),
                                          atom_offsets_.end(), pdb_atom_id);
            pdb_atom_id -= lower - atom_offsets_.begin();
            return static_cast<size_t>(pdb_atom_id) - atom_offsets_.front();
        } catch (const Error&) {
            throw format_error("could not read atomic number in '{}'", line);
        }
    };

    auto i = read_index(6);

    if (line_length > 11) {
        auto j = read_index(11);
        add_bond(i, j);
    } else {
        return;
    }

    if (line_length > 16) {
        auto j = read_index(16);
        add_bond(i, j);
    } else {
        return;
    }

    if (line_length > 21) {
        auto j = read_index(21);
        add_bond(i, j);
    } else {
        return;
    }

    if (line_length > 26) {
        auto j = read_index(26);
        add_bond(i, j);
    } else {
        return;
    }
}

void PDBFormat::chain_ended(Frame& frame) {
    for (const auto& residue: residues_) {
        frame.add_residue(residue.second);
    }

    // This is a 'hack' to allow for badly formatted PDB files which restart
    // the residue ID after a TER residue in cases where they should not.
    // IE a metal Ion given the chain ID of A and residue ID of 1 even though
    // this residue already exists.
    residues_.clear();
}

void PDBFormat::link_standard_residue_bonds(Frame& frame) {
    bool link_previous_peptide = false;
    bool link_previous_nucleic = false;
    int64_t previous_residue_id = 0;
    size_t previous_carboxylic_id = 0;

    for (const auto& residue: frame.topology().residues()) {
        auto residue_table = PDBConnectivity::find(residue.name());
        if (!residue_table) {
            continue;
        }

        std::map<std::string, size_t> atom_name_to_index;
        for (size_t atom : residue) {
            atom_name_to_index[frame[atom].name()] =  atom;
        }

        const auto& amide_nitrogen = atom_name_to_index.find("N");
        const auto& amide_carbon = atom_name_to_index.find("C");

        if (!residue.id()) {
            warning("PDB reader", "got a residues without id, this should not happen");
            continue;
        }

        auto resid = *residue.id();
        if (link_previous_peptide &&
            amide_nitrogen != atom_name_to_index.end() &&
            resid == previous_residue_id + 1 )
        {
            link_previous_peptide = false;
            frame.add_bond(previous_carboxylic_id, amide_nitrogen->second);
        }

        if (amide_carbon != atom_name_to_index.end() ) {
            link_previous_peptide = true;
            previous_carboxylic_id = amide_carbon->second;
            previous_residue_id = resid;
        }

        const auto& three_prime_oxygen = atom_name_to_index.find("O3'");
        const auto& five_prime_phosphorus = atom_name_to_index.find("P");

        if (link_previous_nucleic &&
            five_prime_phosphorus != atom_name_to_index.end() &&
            resid == previous_residue_id + 1 )
        {
            link_previous_nucleic = false;
            frame.add_bond(previous_carboxylic_id, three_prime_oxygen->second);
        }

        if (three_prime_oxygen != atom_name_to_index.end() ) {
            link_previous_nucleic = true;
            previous_carboxylic_id = three_prime_oxygen->second;
            previous_residue_id = resid;
        }

        // A special case missed by the standards committee????
        if (atom_name_to_index.count("HO5'") != 0) {
            frame.add_bond(atom_name_to_index["HO5'"], atom_name_to_index["O5'"]);
        }

        for (const auto& link: *residue_table) {
            const auto& first_atom = atom_name_to_index.find(link.first);
            const auto& second_atom = atom_name_to_index.find(link.second);

            if (first_atom == atom_name_to_index.end()) {
                const auto& first_name = link.first.string();
                if (first_name[0] != 'H' && first_name != "OXT" &&
                    first_name[0] != 'P' && first_name.substr(0, 2) != "OP" ) {
                    warning("PDB reader",
                        "could not find standard atom '{}' in residue '{}' (resid {})",
                        first_name, residue.name(), resid
                    );
                }
                continue;
            }

            if (second_atom == atom_name_to_index.end()) {
                const auto& second_name = link.second.string();
                if (second_name[0] != 'H' && second_name != "OXT" &&
                    second_name[0] != 'P' && second_name.substr(0, 2) != "OP" ) {
                        warning("PDB reader",
                            "could not find standard atom '{}' in residue '{}' (resid {})",
                            second_name, residue.name(), resid
                        );
                }
                continue;
            }

            frame.add_bond(first_atom->second, second_atom->second);
        }
    }
}

Record get_record(std::string_view line) {
    auto rec = line.substr(0, 6);
    if(rec == "ENDMDL") {
        return Record::ENDMDL;
    } else if(rec.substr(0, 3) == "END") {
        // Handle missing whitespace in END record
        return Record::END;
    } else if (rec == "CRYST1") {
        return Record::CRYST1;
    } else if (rec == "ATOM  ") {
        return Record::ATOM;
    } else if (rec == "HETATM") {
        return Record::HETATM;
    } else if (rec == "CONECT") {
        return Record::CONECT;
    } else if (rec.substr(0, 5) == "MODEL") {
        return Record::MODEL;
    } else if (rec.substr(0, 3) == "TER") {
        return Record::TER;
    } else if (rec == "HELIX ") {
        return Record::HELIX;
    } else if (rec == "SHEET ") {
        return Record::SHEET;
    } else if (rec == "TURN  ") {
        return Record::TURN;
    } else if (rec == "HEADER") { // These appear the least, so check last
        return Record::HEADER;
    } else if (rec == "TITLE ") {
        return Record::TITLE;
    } else if (rec == "REMARK" || rec == "MASTER" || rec == "AUTHOR" ||
               rec == "CAVEAT" || rec == "COMPND" || rec == "EXPDTA" ||
               rec == "KEYWDS" || rec == "OBSLTE" || rec == "SOURCE" ||
               rec == "SPLIT " || rec == "SPRSDE" || rec == "JRNL  " ||
               rec == "SEQRES" || rec == "HET   " || rec == "REVDAT" ||
               rec == "SCALE1" || rec == "SCALE2" || rec == "SCALE3" ||
               rec == "ORIGX1" || rec == "ORIGX2" || rec == "ORIGX3" ||
               rec == "SCALE1" || rec == "SCALE2" || rec == "SCALE3" ||
               rec == "ANISOU" || rec == "SITE  " || rec == "FORMUL" ||
               rec == "DBREF " || rec == "HETNAM" || rec == "HETSYN" ||
               rec == "SSBOND" || rec == "LINK  " || rec == "SEQADV" ||
               rec == "MODRES" || rec == "SEQRES" || rec == "CISPEP") {
        return Record::IGNORED_;
    } else if (trim(line).empty()) {
        return Record::IGNORED_;
    } else {
        return Record::UNKNOWN_;
    }
}

static std::string to_pdb_index(int64_t value, size_t width) {
    auto encoded = encode_hybrid36(width, value + 1);

    if (encoded[0] == '*' && (value == MAX_HYBRID36_W4_NUMBER || value == MAX_HYBRID36_W5_NUMBER)) {
        const auto* type = width == 5 ? "atom" : "residue";
        warning("PDB writer", "the value for a {} serial/id is too large, using '{}' instead", type, encoded);
    }

    return encoded;
}

struct ResidueInformation {
    std::string atom_hetatm = "HETATM";
    std::string resname;
    std::string resid;
    std::string chainid;
    std::string insertion_code;
    std::string composition_type;
    std::string segment;
};

static ResidueInformation get_residue_information(optional<const Residue&> residue_opt, int64_t& max_resid) {
    ResidueInformation info;

    if (!residue_opt) {
        auto value = max_resid++;
        info.resid = to_pdb_index(value, 4);
        return info;
    }

    const auto& residue = *residue_opt;

    if (residue.get<Property::BOOL>("is_standard_pdb").value_or(false)) {
        // only use ATOM if the residue is standardized
        info.atom_hetatm = "ATOM  ";
    }

    info.resname = residue.name();
    if (info.resname.length() > 3) {
        warning("PDB writer", "residue '{}' name is too long, it will be truncated", info.resname);
        info.resname = info.resname.substr(0, 3);
    }

    if (residue.id()) {
        info.resid = to_pdb_index(residue.id().value() - 1, 4);
    }

    info.chainid = residue.get<Property::STRING>("chainid").value_or(" ");
    if (info.chainid.length() > 1) {
        warning("PDB writer",
            "residue's chain id '{}' is too long, it will be truncated",
            info.chainid
        );
        info.chainid = info.chainid[0];
    }


    info.insertion_code = residue.get<Property::STRING>("insertion_code").value_or("");
    if (info.insertion_code.length() > 1) {
        warning("PDB writer",
            "residue's insertion code '{}' is too long, it will be truncated",
            info.insertion_code
        );
        info.insertion_code = info.insertion_code[0];
    }

    info.segment = residue.get<Property::STRING>("segname").value_or("");
    if (info.segment.length() > 4) {
        warning("PDB writer",
            "residue's segment name '{}' is too long, it will be truncated",
            info.segment
        );
        info.segment = info.segment.substr(0, 4);
    }

    info.composition_type = residue.get<Property::STRING>("composition_type").value_or("");

    return info;
}

static bool needs_ter_record(const ResidueInformation& residue) {
    return !(
        residue.composition_type.empty() ||
        residue.composition_type == "other" ||
        residue.composition_type == "OTHER" ||
        residue.composition_type == "non-polymer" ||
        residue.composition_type == "NON-POLYMER"
    );
}

// This function adjusts a given index to account for intervening TER records.
// It does so by determining the position of the greatest TER record in `ters`
// and uses iterator arithmetic to calculate the adjustment. Note that `ters`
// is expected to be sorted
static int64_t adjust_for_ter_residues(size_t v, const std::vector<size_t>& ters) {
    auto lower = std::lower_bound(ters.begin(), ters.end(), v + 1);
    auto b0 = static_cast<int64_t>(v) + (lower - ters.begin());
    return b0;
}

void PDBFormat::write_next(const Frame& frame) {
    written_ = true;
    file_.print("MODEL {:>4}\n", models_ + 1);

    auto lengths = frame.cell().lengths();
    auto angles = frame.cell().angles();
    check_values_size(lengths, 9, "cell lengths");
    check_values_size(angles, 7, "cell angles");
    // Do not try to guess the space group and the z value, just use the
    // default one.
    file_.print("CRYST1{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f} P 1           1\n",
        lengths[0], lengths[1], lengths[2], angles[0], angles[1], angles[2]
    );

    // Only use numbers bigger than the biggest residue id as "resSeq" for
    // atoms without associated residue.
    int64_t max_resid = 0;
    for (const auto& residue: frame.topology().residues()) {
        auto resid = residue.id();
        if (resid && resid.value() > max_resid) {
            max_resid = resid.value();
        }
    }

    // Used to skip writing unnecessary connect record
    // use std::deque because std::vector<bool> have a surprising behavior due
    // to C++ standard requiring it to pack multiple bool in a byte
    auto is_atom_record = std::deque<bool>(frame.size(), false);

    // Used for writing TER records.
    size_t ter_count = 0;
    optional<ResidueInformation> last_residue = nullopt;
    std::vector<size_t> ter_serial_numbers;

    const auto& positions = frame.positions();
    for (size_t i = 0; i < frame.size(); i++) {

        auto altloc = frame[i].get<Property::STRING>("altloc").value_or(" ");
        if (altloc.length() > 1) {
            warning("PDB writer", "altloc '{}' is too long, it will be truncated", altloc);
            altloc = altloc[0];
        }

        auto residue = frame.topology().residue_for_atom(i);
        auto resinfo = get_residue_information(residue, max_resid);
        if (resinfo.atom_hetatm == "ATOM  ") {
            is_atom_record[i] = true;
        }

        assert(resinfo.resname.length() <= 3);

        if (last_residue && last_residue->chainid != resinfo.chainid && needs_ter_record(*last_residue)) {
            file_.print("TER   {: >5}      {:3} {:1}{: >4s}{:1}\n",
                to_pdb_index(static_cast<int64_t>(i + ter_count), 5),
                last_residue->resname, last_residue->chainid, last_residue->resid, last_residue->insertion_code);
            ter_serial_numbers.push_back(i + ter_count);
            ++ter_count;
        }

        const auto& pos = positions[i];
        check_values_size(pos, 8, "atomic position");
        file_.print(
            "{: <6}{: >5} {: <4s}{:1}{:3} {:1}{: >4s}{:1}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}      {: <4s}{: >2s}\n",
            resinfo.atom_hetatm, to_pdb_index(static_cast<int64_t>(i + ter_count), 5), frame[i].name(), altloc,
            resinfo.resname, resinfo.chainid, resinfo.resid, resinfo.insertion_code,
            pos[0], pos[1], pos[2], 1.0, 0.0, resinfo.segment, frame[i].type()
        );

        if (residue) {
            last_residue = std::move(resinfo);
        } else {
            last_residue = nullopt;
        }

    }

    auto connect = std::vector<std::vector<int64_t>>(frame.size());
    for (const auto& bond : frame.topology().bonds()) {
        if (is_atom_record[bond[0]] && is_atom_record[bond[1]]) { // both must be standard residue atoms
            continue;
        }
        if (bond[0] > 87440031 || bond[1] > 87440031) {
            warning("PDB writer", "atomic index is too big for CONNECT, removing the bond between {} and {}", bond[0], bond[1]);
            continue;
        }

        connect[bond[0]].push_back(adjust_for_ter_residues(bond[1], ter_serial_numbers));
        connect[bond[1]].push_back(adjust_for_ter_residues(bond[0], ter_serial_numbers));
    }

    for (size_t i = 0; i < frame.size(); i++) {
        auto connections = connect[i].size();
        auto lines = connections / 4 + 1;
        if (connections == 0) {
            continue;
        }

        auto correction = adjust_for_ter_residues(i, ter_serial_numbers);

        for (size_t conect_line = 0; conect_line < lines; conect_line++) {
            file_.print("CONECT{: >5}", to_pdb_index(correction, 5));

            auto last = std::min(connections, 4 * (conect_line + 1));
            for (size_t j = 4 * conect_line; j < last; j++) {
                file_.print("{: >5}", to_pdb_index(connect[i][j], 5));
            }
            file_.print("\n");
        }
    }

    file_.print("ENDMDL\n");

    models_++;
}

void check_values_size(const Vector3D& values, unsigned width, const std::string& context) {
    double max_pos = std::pow(10.0, width) - 1;
    double max_neg = -std::pow(10.0, width - 1) + 1;
    if (values[0] > max_pos || values[1] > max_pos || values[2] > max_pos ||
        values[0] < max_neg || values[1] < max_neg || values[2] < max_neg) {
        throw format_error(
            "value in {} is too big for representation in PDB format", context
        );
    }
}

PDBFormat::~PDBFormat() {
    if (written_) {
      file_.print("END\n");
    }
}

optional<uint64_t> PDBFormat::forward() {
    auto position = file_.tellpos();

    while (!file_.eof()) {
        auto line = file_.readline();

        std::string saved_line;
        if (line.substr(0, 6) == "ENDMDL") {
            // save the current line into an owned string, since the underlying
            // buffer in file_ can change in file_->readline();
            saved_line = std::string(line);
            line = saved_line;

            auto save = file_.tellpos();
            auto next = file_.readline();
            file_.seekpos(save);

            if (next.substr(0, 3) == "END") {
                // We found another record starting by END in the next line,
                // we skip this one and wait for the next one
                continue;
            }
        }

        if (line.substr(0, 3) == "END") {
            return position;
        }
    }

    // Handle file without END/ENDML record at all
    if (position == 0) {
        return position;
    } else {
        return nullopt;
    }
}
