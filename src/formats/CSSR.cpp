// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cstddef>
#include <cstdint>

#include <array>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "chemfiles/types.hpp"
#include "chemfiles/parse.hpp"
#include "chemfiles/utils.hpp"
#include "chemfiles/warnings.hpp"
#include "chemfiles/error_fmt.hpp"
#include "chemfiles/external/optional.hpp"

#include "chemfiles/File.hpp"
#include "chemfiles/Format.hpp"
#include "chemfiles/Atom.hpp"
#include "chemfiles/Frame.hpp"
#include "chemfiles/Topology.hpp"
#include "chemfiles/UnitCell.hpp"
#include "chemfiles/Connectivity.hpp"
#include "chemfiles/FormatMetadata.hpp"

#include "chemfiles/formats/CSSR.hpp"

namespace chemfiles {
    class MemoryBuffer;
}

using namespace chemfiles;

template<> const FormatMetadata& chemfiles::format_metadata<CSSRFormat>() {
    static FormatMetadata metadata;
    metadata.name = "CSSR";
    metadata.extension = ".cssr";
    metadata.description = "CSSR text format";
    metadata.reference = "http://www.chem.cmu.edu/courses/09-560/docs/msi/modenv/D_Files.html#944777";

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


CSSRFormat::CSSRFormat(std::string path, File::Mode mode, File::Compression compression)
    : TextFormat(std::move(path), mode, compression)
{
    if (mode == File::APPEND) {
        throw format_error("append mode ('a') is not supported with CSSR format");
    }
}

CSSRFormat::CSSRFormat(std::shared_ptr<MemoryBuffer> memory, File::Mode mode, File::Compression compression)
    : TextFormat(std::move(memory), mode, compression)
{
    if (mode == File::APPEND) {
        throw format_error("append mode ('a') is not supported with CSSR format");
    }
}

void CSSRFormat::read_next(Frame& frame) {
    if (file_.tellpos() != 0) {
        throw format_error("CSSR format only supports reading one frame");
    }

    // Read unit cell
    Vector3D lengths;
    Vector3D angles;
    scan(file_.readline().substr(38), lengths[0], lengths[1], lengths[2]);
    scan(file_.readline().substr(21), angles[0], angles[1], angles[2]);
    frame.set_cell({lengths, angles});

    size_t natoms = 0;
    int coordinate_style = -1;
    scan(file_.readline(), natoms, coordinate_style);
    bool use_fractional = (coordinate_style == 0);

    // Title line
    file_.readline();

    frame.reserve(natoms);
    std::vector<std::vector<size_t>> connectivity(natoms);
    for (size_t i=0; i<natoms; i++) {
        auto line = file_.readline();
        unsigned atom_id = 0;
        std::string name;
        double x = 0;
        double y = 0;
        double z = 0;
        unsigned bonds[8] = {0};
        double charge = 0;

        scan(line,
            atom_id, name, x, y, z,
            bonds[0], bonds[1], bonds[2], bonds[3],
            bonds[4], bonds[5], bonds[6], bonds[7],
            charge
        );

        auto position = Vector3D(x, y, z);
        if (use_fractional) {
            position = frame.cell().matrix() * position;
        }

        // Atomic names can be created as <type><id>: O121 H22
        auto type = name;
        size_t type_length = 0;
        for (auto ch: type) {
            if (is_ascii_digit(ch)) {
                break;
            }
            type_length += 1;
        }
        type = type.substr(0, type_length);

        auto atom = Atom(std::move(name), std::move(type));
        atom.set_charge(charge);
        frame.add_atom(std::move(atom), position);

        for (auto bond: bonds) {
            if (bond != 0) {
                connectivity[atom_id - 1].push_back(bond - 1);
            }
        }
    }

    for (size_t i=0; i<natoms; i++) {
        for (auto j: connectivity[i]) {
            frame.add_bond(i, j);
        }
    }
}

void CSSRFormat::write_next(const Frame& frame) {
    if (file_.tellpos() != 0) {
        throw format_error("CSSR format only supports writing one frame");
    }

    auto lengths = frame.cell().lengths();
    file_.print(
        " REFERENCE STRUCTURE = 00000   A,B,C ={:8.3f}{:8.3f}{:8.3f}\n",
        lengths[0], lengths[1], lengths[2]
    );

    auto angles = frame.cell().angles();
    file_.print(
        "   ALPHA,BETA,GAMMA ={:8.3f}{:8.3f}{:8.3f}    SPGR =  1 P1\n",
        angles[0], angles[1], angles[2]
    );

    if (frame.size() > 9999) {
        warning("CCSR writer", "too many atoms, the file might not open with other programs");
        file_.print("{} 0\n", frame.size());
    } else {
        file_.print("{:4}   0\n", frame.size());
    }

    // TODO: use the frame name/title property in the file title
    file_.print(" file created with chemfiles\n", frame.size());

    auto connectivity = std::vector<std::vector<size_t>>(frame.size());
    for (const auto& bond : frame.topology().bonds()) {
        if (bond[0] > 9999 || bond[1] > 9999) {
            warning("CCSR writer", "atomic index is too big for connectivity record, removing the bond");
            continue;
        }
        connectivity[bond[0]].push_back(bond[1]);
        connectivity[bond[1]].push_back(bond[0]);
    }

    const auto& positions = frame.positions();
    auto cell_inv = frame.cell().matrix().invert();
    for (size_t i = 0; i<frame.size(); i++) {
        std::string atom_id;
        if (i <= 9999) {
            atom_id = std::to_string(i + 1);
        } else {
            atom_id = "****";
        }

        auto fractional = cell_inv * positions[i];
        file_.print(
            "{:4} {:4}  {:9.5f} {:9.5f} {:9.5f}",
            atom_id, frame[i].name(), fractional[0], fractional[1], fractional[2]
        );

        auto bonds = 0;
        for (auto bond: connectivity[i]) {
            if (bonds >= 8) {
                warning("CCSR writer", "too many bonds with atom {}, only 8 are supported", i);
                break;
            }
            file_.print("{:4}", bond + 1);
            bonds += 1;
        }
        while (bonds < 8) {
            file_.print("   0");
            bonds += 1;
        }

        file_.print(" {:7.3f}\n", frame[i].charge());
    }
}

optional<uint64_t> CSSRFormat::forward() {
    // CSSR only supports one step, so always act like there is only one
    auto position = file_.tellpos();
    if (position == 0) {
        // advance the pointer so that the next call to forward returns -1
        file_.readline();
        return position;
    } else {
        return nullopt;
    }
}
