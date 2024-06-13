// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license
#include "catch.hpp"
#include "chemfiles/Connectivity.hpp"
#include "helpers.hpp"
#include "chemfiles.hpp"
using namespace chemfiles;

#include "fmt/core.h"


TEST_CASE("Read files in LAMMPS data format") {
    SECTION("File created by VMD/Topotools") {
        auto file = Trajectory("data/lammps-data/solvated.lmp", 'r', "LAMMPS Data");
        REQUIRE(file.nsteps() == 1);

        auto frame = file.read();
        CHECK(frame.size() == 7772);

        CHECK(frame.cell().shape() == UnitCell::ORTHORHOMBIC);
        CHECK(approx_eq(frame.cell().lengths(), {34.023998, 34.023998, 163.035995}, 1e-12));

        auto positions = frame.positions();
        CHECK(approx_eq(positions[0], {4.253000, 12.759000, 63.506001}, 1e-12));
        CHECK(approx_eq(positions[364], {8.134000, 2.322000, 82.219002}, 1e-12));
        CHECK(approx_eq(positions[653], {6.184000, 8.134000, 104.334000}, 1e-12));

        auto& topology = frame.topology();
        CHECK(topology.bonds().size() == 6248);
        CHECK(topology[3].name() == "Zn");
        CHECK(topology[3].type() == "Zn");
        CHECK(topology[3].mass() == 65.408997);

        CHECK(topology[12].name() == "C1");
        CHECK(topology[12].type() == "C1");
        CHECK(topology[12].mass() == 42.0);

        // Check the read_step function
        frame = file.read_step(0);
        positions = frame.positions();
        CHECK(approx_eq(positions[0], {4.253000, 12.759000, 63.506001}, 1e-12));
        CHECK(approx_eq(positions[364], {8.134000, 2.322000, 82.219002}, 1e-12));
        CHECK(approx_eq(positions[653], {6.184000, 8.134000, 104.334000}, 1e-12));
    }

    SECTION("File created with LAMMPS") {
        auto file = Trajectory("data/lammps-data/data.body", 'r', "LAMMPS Data");
        Frame frame = file.read();

        CHECK(frame.size() == 100);
        CHECK(frame.cell().shape() == UnitCell::ORTHORHOMBIC);
        CHECK(approx_eq(frame.cell().lengths(), {31.064449134, 31.064449134, 1.0}, 1e-12));

        auto positions = frame.positions();
        CHECK(approx_eq(positions[0], {-15.5322, -15.5322, 0.0}, 1e-12));
        CHECK(approx_eq(positions[22], {-9.31933, -9.31933, 0.0}, 1e-12));

        REQUIRE(frame.velocities());
        auto velocities = *frame.velocities();
        CHECK(approx_eq(velocities[0], {1.02255489961, 2.92322463726, 4.88805110017}));
        CHECK(approx_eq(velocities[1], {0.111646059519, 0.474226666855, 0.68604865644}));
        CHECK(approx_eq(velocities[5], {1.14438145745, 4.42784814304, 1.75516442452}));
        CHECK(approx_eq(velocities[42], {4.70147770939, 2.13317266836, 1.29333445263}));

        auto& topology = frame.topology();
        CHECK(topology.bonds().size() == 0);
        CHECK(topology[0].mass() == 6);
        CHECK(topology[1].mass() == 4);
        CHECK(topology[2].mass() == 3);

        CHECK(topology[12].name() == "1");
        CHECK(topology[12].type() == "1");
    }

    SECTION("Triclinic cells") {
        auto file = Trajectory("data/lammps-data/triclinic-1.lmp", 'r', "LAMMPS Data");
        auto cell = file.read().cell();

        CHECK(cell.shape() == UnitCell::TRICLINIC);
        CHECK(cell.lengths() == Vector3D(34, 34, 34));
        CHECK(cell.angles() == Vector3D(90, 90, 90));

        file = Trajectory("data/lammps-data/triclinic-2.lmp", 'r', "LAMMPS Data");
        cell = file.read().cell();
        CHECK(cell.shape() == UnitCell::TRICLINIC);
        CHECK(approx_eq(cell.lengths(), {34, 34.3656805549, 35.0570962859}, 1e-9));
        CHECK(approx_eq(cell.angles(), {87.0501134427, 103.1910720469, 81.634113876}, 1e-9));
    }

    SECTION("Molecule ids") {
        auto file = Trajectory("data/lammps-data/molid.lmp", 'r', "LAMMPS Data");
        auto frame = file.read();

        CHECK(frame.size() == 12);
        CHECK(frame.topology().residues().size() == 3);

        auto& topology = frame.topology();
        CHECK_FALSE(topology.residue_for_atom(0));
        CHECK_FALSE(topology.residue_for_atom(1));
        CHECK_FALSE(topology.residue_for_atom(2));

        CHECK(topology.residue_for_atom(3)->contains(4));
        CHECK(topology.residue_for_atom(3)->contains(5));
        CHECK(topology.residue_for_atom(3)->id().value() == 1);

        CHECK(topology.residue_for_atom(6)->contains(7));
        CHECK(topology.residue_for_atom(6)->contains(8));
        CHECK(topology.residue_for_atom(6)->id().value() == 2);

        CHECK(topology.residue_for_atom(9)->contains(10));
        CHECK(topology.residue_for_atom(9)->contains(11));
        CHECK(topology.residue_for_atom(9)->id().value() == 3);
    }

    SECTION("File with whitespaces") {
        // https://github.com/chemfiles/chemfiles/issues/485
        auto file = Trajectory("data/lammps-data/whitespaces.lmp", 'r', "LAMMPS Data");
        auto frame = file.read();
    }
}

TEST_CASE("Write files in LAMMPS data format") {
    const auto EXPECTED_CONTENT =
    "LAMMPS data file -- atom_style full -- generated by chemfiles\n\n"
    "6 atoms\n"
    "4 bonds\n"
    "4 angles\n"
    "2 dihedrals\n"
    "1 impropers\n"
    "4 atom types\n"
    "3 bond types\n"
    "3 angle types\n"
    "2 dihedral types\n"
    "1 improper types\n"
    "0.00000000 5.00000000 xlo xhi\n"
    "0.00000000 6.06217783 ylo yhi\n"
    "0.00000000 9.00000000 zlo zhi\n"
    "1.50000000 0.00000000 0.00000000 xy xz yz\n"
    "\n"
    "# Pair Coeffs\n"
    "# 1 As\n"
    "# 2 As\n"
    "# 3 B\n"
    "# 4 C\n"
    "\n"
    "# Bond Coeffs\n"
    "# 1 As-B\n"
    "# 2 B-B\n"
    "# 3 B-C\n"
    "\n"
    "# Angle Coeffs\n"
    "# 1 As-B-B\n"
    "# 2 As-B-C\n"
    "# 3 B-B-C\n"
    "\n"
    "# Dihedrals Coeffs\n"
    "# 1 As-B-B-C\n"
    "# 2 C-B-B-C\n"
    "\n"
    "# Impropers Coeffs\n"
    "# 1 As-B-B-C\n"
    "\n"
    "Masses\n"
    "\n"
    "1 25.0000 # As\n"
    "2 74.9216 # As\n"
    "3 10.8100 # B\n"
    "4 12.0110 # C\n"
    "\n"
    "Atoms # full\n"
    "\n"
    "1 1 1 0.00000 1.10000 2.20000 3.30000 # As\n"
    "2 2 2 0.00000 1.10000 2.20000 3.30000 # As\n"
    "3 2 3 -2.40000 1.10000 2.20000 3.30000 # B\n"
    "4 2 4 0.00000 1.10000 2.20000 3.30000 # C\n"
    "5 2 3 0.00000 1.10000 2.20000 3.30000 # B\n"
    "6 2 4 0.00000 1.10000 2.20000 3.30000 # C\n"
    "\n"
    "Velocities\n"
    "\n"
    "1 0.1 0.2 0.3\n"
    "2 0.1 0.2 0.3\n"
    "3 0.1 0.2 0.3\n"
    "4 0.1 0.2 0.3\n"
    "5 0.1 0.2 0.3\n"
    "6 0.1 0.2 0.3\n"
    "\n"
    "Bonds\n"
    "\n"
    "1 1 2 3\n"
    "2 3 3 4\n"
    "3 2 3 5\n"
    "4 3 5 6\n"
    "\n"
    "Angles\n"
    "\n"
    "1 2 2 3 4\n"
    "2 1 2 3 5\n"
    "3 3 3 5 6\n"
    "4 3 4 3 5\n"
    "\n"
    "Dihedrals\n"
    "\n"
    "1 1 2 3 5 6\n"
    "2 2 4 3 5 6\n"
    "\n"
    "Impropers\n"
    "\n"
    "1 1 2 3 4 5\n";

    auto tmpfile = NamedTempPath(".lmp");

    auto frame = Frame(UnitCell({5, 7, 9}, {90, 90, 120}));
    frame.add_velocities();
    frame.add_atom(Atom("As"), {1.1, 2.2, 3.3}, {0.1, 0.2, 0.3});
    frame.add_atom(Atom("As"), {1.1, 2.2, 3.3}, {0.1, 0.2, 0.3});
    frame.add_atom(Atom("B"), {1.1, 2.2, 3.3}, {0.1, 0.2, 0.3});
    frame.add_atom(Atom("C"), {1.1, 2.2, 3.3}, {0.1, 0.2, 0.3});
    frame.add_atom(Atom("B"), {1.1, 2.2, 3.3}, {0.1, 0.2, 0.3});
    frame.add_atom(Atom("C"), {1.1, 2.2, 3.3}, {0.1, 0.2, 0.3});
    frame.add_bond(2, 1);
    frame.add_bond(2, 3);
    frame.add_bond(2, 4);
    frame.add_bond(4, 5);

    frame[0].set_mass(25);
    frame[2].set_charge(-2.4);

    auto trajectory = Trajectory(tmpfile, 'w', "LAMMPS Data");
    trajectory.write(frame);

    CHECK_THROWS_WITH(
        trajectory.write(frame),
        "LAMMPS Data format only supports writing one frame"
    );
    trajectory.close();

    auto content = read_text_file(tmpfile);
    CHECK(content == EXPECTED_CONTENT);
}

TEST_CASE("Write files in LAMMPS data format with explicit types") {
    const auto EXPECTED_CONTENT =
    "LAMMPS data file -- atom_style full -- generated by chemfiles\n\n"
    "6 atoms\n"
    "4 bonds\n"
    "4 angles\n"
    "2 dihedrals\n"
    "1 impropers\n"
    "4 atom types\n"
    "4 bond types\n"
    "4 angle types\n"
    "2 dihedral types\n"
    "1 improper types\n"
    "0.00000000 5.00000000 xlo xhi\n"
    "0.00000000 6.06217783 ylo yhi\n"
    "0.00000000 9.00000000 zlo zhi\n"
    "1.50000000 0.00000000 0.00000000 xy xz yz\n"
    "\n"
    "# Pair Coeffs\n"
    "# 1 As\n"
    "# 2 As\n"
    "# 3 B\n"
    "# 4 C\n"
    "\n"
    "# Bond Coeffs\n"
    "# 1 As-B\n"
    "# 2 B-B\n"
    "# 3 B-C\n"
    "# 4 B-C\n"
    "\n"
    "# Angle Coeffs\n"
    "# 1 As-B-B\n"
    "# 2 As-B-C\n"
    "# 3 B-B-C\n"
    "# 4 B-B-C\n"
    "\n"
    "# Dihedrals Coeffs\n"
    "# 1 As-B-B-C\n"
    "# 2 C-B-B-C\n"
    "\n"
    "# Impropers Coeffs\n"
    "# 1 As-B-B-C\n"
    "\n"
    "Masses\n"
    "\n"
    "1 25.0 # As\n"
    "2 74.9216 # As\n"
    "3 10.81 # B\n"
    "4 12.011 # C\n"
    "\n"
    "Atoms # full\n"
    "\n"
    "1 1 1 0.0 1.1 2.2 3.3 # As\n"
    "2 2 2 0.0 1.1 2.2 3.3 # As\n"
    "3 2 3 -2.4 1.1 2.2 3.3 # B\n"
    "4 2 4 0.0 1.1 2.2 3.3 # C\n"
    "5 2 3 0.0 1.1 2.2 3.3 # B\n"
    "6 2 4 0.0 1.1 2.2 3.3 # C\n"
    "\n"
    "Velocities\n"
    "\n"
    "1 0.1 0.2 0.3\n"
    "2 0.1 0.2 0.3\n"
    "3 0.1 0.2 0.3\n"
    "4 0.1 0.2 0.3\n"
    "5 0.1 0.2 0.3\n"
    "6 0.1 0.2 0.3\n"
    "\n"
    "Bonds\n"
    "\n"
    "1 1 2 3\n"
    "2 2 3 4\n"
    "3 3 3 5\n"
    "4 4 5 6\n"
    "\n"
    "Angles\n"
    "\n"
    "1 1 2 3 4\n"
    "2 2 2 3 5\n"
    "3 3 3 5 6\n"
    "4 4 4 3 5\n"
    "\n"
    "Dihedrals\n"
    "\n"
    "1 1 2 3 5 6\n"
    "2 2 4 3 5 6\n"
    "\n"
    "Impropers\n"
    "\n"
    "1 1 2 3 4 5\n";

    auto tmpfile = NamedTempPath(".lmp");

    auto frame = Frame(UnitCell({5, 7, 9}, {90, 90, 120}));
    frame.add_velocities();
    frame.add_atom(Atom("As"), {1.1, 2.2, 3.3}, {0.1, 0.2, 0.3});
    frame.add_atom(Atom("As"), {1.1, 2.2, 3.3}, {0.1, 0.2, 0.3});
    frame.add_atom(Atom("B"), {1.1, 2.2, 3.3}, {0.1, 0.2, 0.3});
    frame.add_atom(Atom("C"), {1.1, 2.2, 3.3}, {0.1, 0.2, 0.3});
    frame.add_atom(Atom("B"), {1.1, 2.2, 3.3}, {0.1, 0.2, 0.3});
    frame.add_atom(Atom("C"), {1.1, 2.2, 3.3}, {0.1, 0.2, 0.3});
    frame.add_bond(2, 1, chemfiles::Bond::UNKNOWN, "1");
    frame.add_bond(2, 3, chemfiles::Bond::UNKNOWN, "2");
    frame.add_bond(2, 4, chemfiles::Bond::UNKNOWN, "3");
    frame.add_bond(4, 5, chemfiles::Bond::UNKNOWN, "4");

    frame.add_angle(1, 2, 3, "1");
    frame.add_angle(1, 2, 4, "2");
    frame.add_angle(2, 4, 5, "3");
    frame.add_angle(3, 2, 4, "4");

    frame.add_dihedral(1, 2, 4, 5, "1");
    frame.add_dihedral(3, 2, 4, 5, "2");

    frame.add_improper(1, 2, 3, 4, "1");

    frame[0].set_mass(25);
    frame[2].set_charge(-2.4);

    auto trajectory = Trajectory(tmpfile, 'w', "LAMMPS Data");
    trajectory.write(frame);

    CHECK_THROWS_WITH(
        trajectory.write(frame),
        "LAMMPS Data format only supports writing one frame"
    );
    trajectory.close();

    auto content = read_text_file(tmpfile);
    CHECK(content == EXPECTED_CONTENT);
}

TEST_CASE("Read and write files in memory") {
    SECTION("Reading from memory") {
        auto content = read_text_file("data/lammps-data/data.body");
        auto file = Trajectory::memory_reader(content.data(), content.size(), "LAMMPS Data");
        auto frame = file.read();

        CHECK(frame.size() == 100);
        CHECK(frame.cell().shape() == UnitCell::ORTHORHOMBIC);
        CHECK(approx_eq(frame.cell().lengths(), {31.064449134, 31.064449134, 1.0}, 1e-12));

        auto positions = frame.positions();
        CHECK(approx_eq(positions[0], {-15.5322, -15.5322, 0.0}, 1e-12));
        CHECK(approx_eq(positions[22], {-9.31933, -9.31933, 0.0}, 1e-12));
    }
}
