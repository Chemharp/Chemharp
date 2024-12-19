// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cassert>
#include <cstddef>
#include <cstdint>

#include <array>
#include <string>
#include <vector>
#include <string_view>

#include "chemfiles/types.hpp"
#include "chemfiles/utils.hpp"
#include "chemfiles/parse.hpp"
#include "chemfiles/error_fmt.hpp"
#include "chemfiles/sorted_set.hpp"
#include "chemfiles/external/optional.hpp"

#include "chemfiles/File.hpp"
#include "chemfiles/Atom.hpp"
#include "chemfiles/Frame.hpp"
#include "chemfiles/Topology.hpp"
#include "chemfiles/UnitCell.hpp"
#include "chemfiles/Connectivity.hpp"
#include "chemfiles/FormatMetadata.hpp"

#include "chemfiles/formats/Tinker.hpp"

using namespace chemfiles;

template<> const FormatMetadata& chemfiles::format_metadata<TinkerFormat>() {
    static FormatMetadata metadata;
    metadata.name = "Tinker";
    metadata.extension = ".arc";
    metadata.description = "Tinker XYZ text format";
    metadata.reference = "http://chembytes.wikidot.com/tnk-tut00#toc2";

    metadata.read = true;
    metadata.write = true;
    metadata.memory = true;

    metadata.positions = true;
    metadata.velocities = false;
    metadata.unit_cell = true;
    metadata.atoms = true;
    metadata.bonds = true;
    metadata.residues = false;
    return metadata;
}

static bool is_unit_cell_line(std::string_view line);

void TinkerFormat::read_next(Frame& frame) {
    auto line = file_.readline();
    // only read the number of atoms, ignore any comment
    size_t n_atoms = 0;
    scan(line, n_atoms);

    std::vector<std::vector<size_t>> bonds(n_atoms);
    frame.reserve(n_atoms);

    auto position = file_.tellpos();
    line = file_.readline();
    if (is_unit_cell_line(line)) {
        // Read the cell
        Vector3D lengths;
        Vector3D angles;
        scan(line, lengths[0], lengths[1], lengths[2], angles[0], angles[1], angles[2]);
        frame.set_cell({lengths, angles});
    } else {
        file_.seekpos(position);
    }

    for (size_t i = 0; i < n_atoms; i++) {
        line = file_.readline();
        double x = 0;
        double y = 0;
        double z = 0;
        int id = 0;
        int atom_type = 0;
        std::string name;
        auto count = scan(line, id, name, x, y, z, atom_type);

        frame.add_atom(Atom(name), Vector3D(x, y, z));
        frame[i].set("atom_type", atom_type);
        while (count != line.size()) {
            size_t bonded = 0;
            count += scan(line.substr(count), bonded);
            bonds[i].push_back(bonded - 1);
        }
    }

    for (size_t i = 0; i < n_atoms; i++) {
        for (size_t j: bonds[i]) {
            frame.add_bond(i, j);
        }
    }
}

void TinkerFormat::write_next(const Frame& frame) {
    auto lengths = frame.cell().lengths();
    auto angles = frame.cell().angles();
    file_.print("{} written by the chemfiles library\n", frame.size());
    file_.print("{:#g} {:#g} {:#g} {:#g} {:#g} {:#g}\n",
        lengths[0], lengths[1], lengths[2], angles[0], angles[1], angles[2]
    );

    const auto& topology = frame.topology();
    // Build type index numbering. type_id will be searched for each atom type,
    // and the position can be used a an unique integer identifier for the atom
    // type.
    auto types_id = sorted_set<std::string>();
    for (const auto& atom: topology) {
        types_id.insert(atom.type());
    }

    // Build bonds index. It will contains all atoms bonded to the atom i in
    // bonded_to[i].
    auto bonded_to = std::vector<std::vector<size_t>>(frame.size());
    for (const auto& bond: topology.bonds()) {
        bonded_to[bond[0]].push_back(bond[1]);
        bonded_to[bond[1]].push_back(bond[0]);
    }

    const auto& positions = frame.positions();
    for (size_t i = 0; i < frame.size(); i++) {
        auto name = topology[i].name();
        if (name.empty()) {
            name = "X";
        }
        auto it = types_id.find(topology[i].type());
        assert(it != types_id.end());
        auto type = (it - types_id.begin()) + 1;

        file_.print("{} {} {:#g} {:#g} {:#g} {}",
            i + 1, name, positions[i][0], positions[i][1], positions[i][2], type
        );
        for (size_t other: bonded_to[i]) {
            file_.print(" {}", other + 1);
        }
        file_.print("\n");
    }
}

optional<uint64_t> TinkerFormat::forward() {
    auto position = file_.tellpos();

    auto line = file_.readline();
    if (trim(line).empty() || file_.eof()) {
        // We just read an empty line, give up here
        return nullopt;
    }

    size_t n_atoms = 0;
    try {
        // Get the number of atoms in the line, ignoring whatever comes after it
        scan(line, n_atoms);
    } catch (const Error&) {
        throw format_error(
            "failed to read number of atoms in '{}' for Tinker XYZ format",
            line
        );
    }

    line = file_.readline();
    // Minus one because we just read a line.
    size_t lines_to_skip = n_atoms - 1;

    if (is_unit_cell_line(line)) {
        lines_to_skip += 1;
    }

    for (size_t i=0; i<lines_to_skip; i++) {
        if (file_.eof()) {
            throw format_error(
                "not enough lines in '{}' for Tinker XYZ format", file_.path()
            );
        }
        file_.readline();
    }

    return position;
}

// This is how tinker does it to check if there is unit cell information
// in the file, so let's follow them here.
bool is_unit_cell_line(std::string_view line) {
    static const char* LETTERS = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
    return line.find_first_of(LETTERS) == std::string::npos;
}
