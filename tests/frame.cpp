// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <catch.hpp>
#include "chemfiles.hpp"
using namespace chemfiles;

TEST_CASE("Frame size", "[Frame]") {
    auto frame = Frame(10);
    CHECK(frame.natoms() == 10);
    CHECK(frame.positions().size() == 10);
    // No velocity data yet
    CHECK_FALSE(frame.velocities());

    frame.resize(15);
    CHECK(frame.natoms() == 15);
    CHECK(frame.positions().size() == 15);
    CHECK_FALSE(frame.velocities());

    frame.add_velocities();
    REQUIRE(frame.velocities());
    CHECK(frame.velocities()->size() == 15);

    frame.resize(30);
    CHECK(frame.natoms() == 30);
    CHECK(frame.positions().size() == 30);
    REQUIRE(frame.velocities());
    CHECK(frame.velocities()->size() == 30);

    frame.resize(2);
    CHECK(frame.natoms() == 2);
    CHECK(frame.positions().size() == 2);
    REQUIRE(frame.velocities());
    CHECK(frame.velocities()->size() == 2);

    frame.add_atom(Atom("H"), vector3d(1, 2, 3), vector3d(4, 5, 6));
    CHECK(frame.natoms() == 3);
    CHECK(frame.positions().size() == 3);
    CHECK(frame.positions()[2] == vector3d(1, 2, 3));
    REQUIRE(frame.velocities());
    CHECK(frame.velocities()->size() == 3);
    CHECK((*frame.velocities())[2] == vector3d(4, 5, 6));

    frame.remove(0);
    CHECK(frame.natoms() == 2);
    CHECK(frame.positions().size() == 2);
    CHECK(frame.velocities()->size() == 2);
}

TEST_CASE("Positions and velocities", "[Frame]") {
    auto frame = Frame(15);

    frame.positions()[0] = vector3d(1, 2, 3);
    CHECK(frame.positions()[0] == vector3d(1, 2, 3));

    frame.add_velocities();
    (*frame.velocities())[0] = vector3d(5, 6, 7);
    CHECK((*frame.velocities())[0] == vector3d(5, 6, 7));

    {
        auto positions = frame.positions();
        auto velocities = frame.velocities();
        for (size_t i=0; i<15; i++) {
            positions[i] = vector3d(4.0, 3.4, 1.0);
            (*velocities)[i] = vector3d(4.0, 3.4, 1.0);
        }
    }

    auto positions = frame.positions();
    auto velocities = frame.velocities();
    for (size_t i=0; i<10; i++){
        CHECK(positions[i] == vector3d(4.0, 3.4, 1.0));
        CHECK((*velocities)[i] == vector3d(4.0, 3.4, 1.0));
    }
}

TEST_CASE("Frame step", "[Frame]") {
    auto frame = Frame();
    CHECK(frame.step() == 0);
    frame.set_step(1000);
    CHECK(frame.step() == 1000);
}

TEST_CASE("Unit cell", "[Frame]") {
    auto frame = Frame();
    CHECK(frame.cell().shape() == UnitCell::INFINITE);
    frame.set_cell(UnitCell(10));
    CHECK(frame.cell().shape() == UnitCell::ORTHORHOMBIC);
}

TEST_CASE("Guess topology", "[Frame]") {
    SECTION("Simple case") {
        auto frame = Frame();
        frame.add_atom(Atom("H"), {{0, 1, 0}});
        frame.add_atom(Atom("O"), {{0, 0, 0}});
        frame.add_atom(Atom("O"), {{1.5, 0, 0}});
        frame.add_atom(Atom("H"), {{1.5, 1, 0}});
        frame.guess_topology();

        auto bonds = std::vector<Bond>{{0, 1}, {1, 2}, {2, 3}};
        auto angles = std::vector<Angle>{{0, 1, 2}, {1, 2, 3}};
        auto dihedrals = std::vector<Dihedral>{{0, 1, 2, 3}};
        CHECK(frame.topology().bonds() == bonds);
        CHECK(frame.topology().angles() == angles);
        CHECK(frame.topology().dihedrals() == dihedrals);
    }

    SECTION("Methane file") {
        auto frame = Trajectory("data/xyz/methane.xyz").read();
        frame.guess_topology();

        auto topology = frame.topology();
        for (size_t i=1; i<5; i++){
            CHECK(topology.isbond(0, i));
        }

        CHECK_FALSE(topology.isbond(2, 4));
        CHECK_FALSE(topology.isbond(1, 2));

        CHECK(topology.isangle(1, 0, 2));
        CHECK(topology.isangle(3, 0, 2));
        CHECK(topology.isangle(2, 0, 4));

        CHECK(topology.bonds().size() == 4);
        CHECK(topology.angles().size() == 6);
        CHECK(topology.dihedrals().size() == 0);

        frame.remove(1);
        topology = frame.topology();
        CHECK(topology.bonds().size() == 3);
        CHECK(topology.angles().size() == 3);

        // Wrong topology size
        frame = Frame(5);
        CHECK_THROWS_AS(frame.set_topology(Topology()), Error);
    }

    SECTION("Cleanup supplementaty H-H bonds") {
        auto frame = Frame();
        frame.add_atom(Atom("O"), {{0, 0, 0}});
        frame.add_atom(Atom("H"), {{0.2, 0.8, 0}});
        frame.add_atom(Atom("H"), {{-0.2, 0.8, 0}});

        frame.guess_topology();
        CHECK(frame.topology().bonds().size() == 2);
        CHECK(frame.topology().isbond(0, 1));
        CHECK(frame.topology().isbond(0, 2));
    }

    // Weird geometries
    SECTION("Triangle molecule") {
        auto frame = Frame();
        frame.add_atom(Atom("C"), {{0, 1, 0}});
        frame.add_atom(Atom("C"), {{0.5, 0, 0}});
        frame.add_atom(Atom("C"), {{-0.5, 0, 0}});

        frame.guess_topology();
        CHECK(frame.topology().bonds().size() == 3);
        CHECK(frame.topology().isbond(0, 1));
        CHECK(frame.topology().isbond(0, 2));
        CHECK(frame.topology().isbond(1, 2));

        CHECK(frame.topology().angles().size() == 3);
        CHECK(frame.topology().isangle(0, 1, 2));
        CHECK(frame.topology().isangle(0, 2, 1));
        CHECK(frame.topology().isangle(1, 0, 2));

        CHECK(frame.topology().dihedrals().size() == 0);
    }

    // Weird geometries
    SECTION("Square molecule") {
        auto frame = Frame();
        frame.add_atom(Atom("C"), {{0, 0, 0}});
        frame.add_atom(Atom("C"), {{1.5, 0, 0}});
        frame.add_atom(Atom("C"), {{1.5, 1.5, 0}});
        frame.add_atom(Atom("C"), {{0, 1.5, 0}});

        frame.guess_topology();
        CHECK(frame.topology().bonds().size() == 4);
        CHECK(frame.topology().isbond(0, 1));
        CHECK(frame.topology().isbond(1, 2));
        CHECK(frame.topology().isbond(2, 3));
        CHECK(frame.topology().isbond(0, 3));

        CHECK(frame.topology().angles().size() == 4);
        CHECK(frame.topology().isangle(0, 1, 2));
        CHECK(frame.topology().isangle(1, 2, 3));
        CHECK(frame.topology().isangle(2, 3, 0));
        CHECK(frame.topology().isangle(3, 0, 1));

        CHECK(frame.topology().dihedrals().size() == 4);
        CHECK(frame.topology().isdihedral(0, 1, 2, 3));
        CHECK(frame.topology().isdihedral(1, 2, 3, 0));
        CHECK(frame.topology().isdihedral(2, 3, 0, 1));
        CHECK(frame.topology().isdihedral(3, 0, 1, 2));
    }
}
