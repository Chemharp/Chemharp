/*
 * Chemharp, an efficient IO library for chemistry file formats
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

#include <sstream>
#include <cassert>

#include "formats/XYZ.hpp"

#include "Error.hpp"
#include "Frame.hpp"
#include "files/File.hpp"

using namespace harp;

std::string XYZFormat::description() const {
    return "XYZ file format.";
}

// Quick forward the file for nsteps
static void forward(TextFile* file, size_t nsteps) {
    size_t i=0;
    // Move the file pointer to the good position
    while (i < nsteps){
        try {
            auto natoms = std::stoul(file->getline());
            file->readlines(natoms+1);
        }
        catch (const FileError& e) {
            throw FormatError("Can not read step " + std::to_string(nsteps) +
                                  ": " + e.what());
        }
        i++;
    }
}

size_t XYZFormat::nsteps(File* file) const {
    auto textfile = dynamic_cast<TextFile*>(file);
    textfile->rewind();
    size_t n = 0;
    while (not textfile->eof()) {
        forward(textfile, 1);
        n++;
    }
    return n;
}

void XYZFormat::read_at_step(File* file, const size_t step, Frame& frame){
    auto textfile = dynamic_cast<TextFile*>(file);
    textfile->rewind();
    forward(textfile, step - 1);
    read_next_step(file, frame);
}

void XYZFormat::read_next_step(File* file, Frame& frame){
    auto textfile = dynamic_cast<TextFile*>(file);
    size_t natoms = std::stoul(textfile->getline());

    textfile->getline(); // XYZ comment line;
    std::vector<std::string> lines;
    lines.reserve(natoms);

    try {
        lines = textfile->readlines(natoms);
    }
    catch (const FileError& e) {
        throw FormatError("Can not read file: " + string(e.what()));
    }

    frame.topology().clear();
    frame.resize(natoms);

    for (size_t i=0; i<lines.size(); i++) {
        std::istringstream string_stream;
        float x, y, z;
        string name;

        string_stream.str(lines[i]);
        string_stream >> name >> x >> y >> z ;
        frame.positions()[i] = Vector3D(x, y, z);
        frame.topology().append(Atom(name));
    }
    frame.topology().guess_bonds();
}

void XYZFormat::write_step(File* file, const Frame& frame){
    auto textfile = dynamic_cast<TextFile*>(file);

    const auto topology = frame.topology();
    const auto positions = frame.positions();
    assert(frame.natoms() == topology.natoms());

    *textfile << frame.natoms() << "\n";
    *textfile << "Written by Chemharp\n";

    for (size_t i=0; i<frame.natoms(); i++){
        auto& pos = positions[i];
        auto name = topology[i].name();
        if (name == "")
            name = "X";
        *textfile << name   << " "
                  << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
    }
}


// Register the xyz format with the ".xyz" extension and the "XYZ" description.
REGISTER(XYZFormat, "XYZ");
REGISTER_EXTENSION(XYZFormat, ".xyz");
