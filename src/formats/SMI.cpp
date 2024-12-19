// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cmath>
#include <cstdlib>
#include <cstdint>

#include <set>
#include <map>
#include <stack>
#include <tuple>
#include <deque>
#include <string>
#include <utility>
#include <vector>
#include <iterator>
#include <algorithm>
#include <string_view>
#include <unordered_map>

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
#include "chemfiles/Connectivity.hpp"
#include "chemfiles/FormatMetadata.hpp"

#include "chemfiles/formats/SMI.hpp"

using namespace chemfiles;

template<> const FormatMetadata& chemfiles::format_metadata<SMIFormat>() {
    static FormatMetadata metadata;
    metadata.name = "SMI";
    metadata.extension = ".smi";
    metadata.description = "SMI text format";
    metadata.reference = "http://opensmiles.org/";

    metadata.read = true;
    metadata.write = true;
    metadata.memory = true;

    metadata.positions = false;
    metadata.velocities = false;
    metadata.unit_cell = false;
    metadata.atoms = true;
    metadata.bonds = true;
    metadata.residues = true;
    return metadata;
}

/// See if all values in array are true
static bool all(const std::deque<bool>& vec);

static bool is_aromatic_organic(char c) {
    switch (c) {
    case 'b':
    case 'c':
    case 'n':
    case 'o':
    case 's':
    case 'p':
        return true;
    default:
        return false;
    }
}

static bool is_aliphatic_organic(std::string_view s) {
    // Note: Including H here is non-standard, but allows for HC(H)(H)H
    static std::set<std::string_view> ALIPHATIC_ORGANIC = {
        "B", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "H"
    };
    return ALIPHATIC_ORGANIC.count(s) != 0;
}

static bool is_chirality_tag(std::string_view s) {
    // list of valid tags
    static std::set<std::string_view> CHIRALITY_TAGS = {
        "TH", "SP", "TB", "OH", "AL"
    };
    return CHIRALITY_TAGS.count(s) != 0;
}

/// Attempts to read a number from the string. The `start` argument is updated
/// to the last numberic character in `smiles`.
static size_t read_number(std::string_view smiles, size_t& start) {
    if (start >= smiles.size() || !is_ascii_digit(smiles[start])) {
        start--;
        return 0;
    }
    size_t old = start;
    do {
        start++;
    } while (start < smiles.size() && is_ascii_digit(smiles[start]));
    start--;
    return parse<size_t>(smiles.substr(old, start - old + 1));
}

static std::string_view read_property_atom(std::string_view smiles, size_t& start) {
    auto old = start;
    if (smiles[start] == '\'') {
        do {
             start++;
        } while( start < smiles.size() && smiles[start] != '\'');
        start++;
        return smiles.substr(old + 1, start - old - 2);
    } else {
        do {
            start++;
        } while (start < smiles.size() && is_ascii_lowercase(smiles[start]));
        return smiles.substr(old, start - old);
    }
}

Atom& SMIFormat::add_atom(Topology& topology, std::string_view atom_name) {
    topology.add_atom(Atom(std::string(atom_name)));

    if (!first_atom_) {
        topology.add_bond(previous_atom_, ++current_atom_, current_bond_order_);
    }

    first_atom_ = false;
    previous_atom_ = current_atom_;
    current_bond_order_ = Bond::SINGLE;
    residues_.back().add_atom(topology.size() - 1);
    return topology[topology.size() - 1];
}

void SMIFormat::process_property_list(Topology& topology, std::string_view smiles) {
    size_t i = 0;
    double mass = 0.0;
    if (is_ascii_digit(smiles[i])) {
        mass = static_cast<double>(read_number(smiles, i));
        ++i;
    }

    bool is_aromatic = is_ascii_lowercase(smiles[i]);
    auto name = read_property_atom(smiles, i);
    auto& new_atom = this->add_atom(topology, name);

    if (is_aromatic) {
        new_atom.set("is_aromatic", true);
    }

    if (mass != 0.0) {
        new_atom.set_mass(mass);
    }

    while (i < smiles.size()) {
        size_t count_or_class = 0;
        std::string chirality_type = "CCW";
        switch (smiles[i]) {
        case 'H':
            count_or_class = read_number(smiles, ++i);
            count_or_class = (count_or_class == 0) ? 1 : count_or_class;
            new_atom.set("hydrogen_count", count_or_class);
            break;
        case '-':
            count_or_class = read_number(smiles, ++i);
            count_or_class = (count_or_class == 0) ? 1 : count_or_class;
            new_atom.set_charge(new_atom.charge() - static_cast<double>(count_or_class));
            break;
        case '+':
            count_or_class = read_number(smiles, ++i);
            count_or_class = (count_or_class == 0) ? 1 : count_or_class;
            new_atom.set_charge(new_atom.charge() + static_cast<double>(count_or_class));
            break;
        case ':':
            count_or_class = read_number(smiles, ++i);
            new_atom.set("smiles_class", count_or_class);
            break;
        case '@':
            // Set the direction of the chirality
            if (i + 1 < smiles.size() && smiles[i + 1] == '@') {
                chirality_type = "CW";
                ++i;
            } else if (i + 2 < smiles.size() && is_chirality_tag(smiles.substr(i + 1, 2))) {
                // Set the type of the chirality
                chirality_type += " " + std::string(smiles.substr(i + 1, 2));
                i += 2;
                chirality_type += std::to_string(read_number(smiles, ++i));
            }
            new_atom.set("chirality", chirality_type);
            break;
        default:
            warning("SMI Reader", "unknown property code: '{}'", smiles[i]);
            break;
        }
        ++i;
    }
}

void SMIFormat::check_ring_(Topology& topology, size_t ring_id) {
    auto ring_lookup = rings_ids_.find(ring_id);

    if (ring_lookup == rings_ids_.end()) {
        rings_ids_.insert({ ring_id, {previous_atom_, current_bond_order_} });
        current_bond_order_ = Bond::SINGLE;
        return;
    }

    // Deviation from the standard, technically bond orders need to be equal,
    // but we will accept the stored order if the current order is single.
    // This is common practice
    topology.add_bond(previous_atom_,
        ring_lookup->second.first,
        current_bond_order_ == Bond::SINGLE ?
        ring_lookup->second.second :
        current_bond_order_
    );
    rings_ids_.erase(ring_lookup);

    current_bond_order_ = Bond::SINGLE;
}

void SMIFormat::read_next(Frame& frame) {
    // Initialize all the reading variables
    branch_point_ = std::stack<size_t, std::vector<size_t>>();
    rings_ids_.clear();
    residues_.clear();
    current_atom_ = 0;
    previous_atom_ = 0;
    current_bond_order_ = Bond::SINGLE;
    first_atom_ = true;

    Topology topology;
    size_t groupid = 1;
    residues_.push_back(Residue("group " + std::to_string(groupid)));

    auto line = file_.readline();
    auto smiles = trim(line);
    while (smiles.empty()) {
        line = file_.readline();
        smiles = trim(line);
    }

    size_t i;
    for (i = 0; i < smiles.size(); ++i) {

        if (is_ascii_whitespace(smiles[i])) {
            break;
        }

        if (smiles[i] == '[') {
            auto prop_end = smiles.find_first_of(']', i);

            if (prop_end == std::string::npos) {
                throw format_error("SMI Reader: unmatched square brace");
            }

            auto square_prop = smiles.substr(i + 1, prop_end - i - 1);
            this->process_property_list(topology, square_prop);

            i = prop_end;

            continue;
        }

        // Check for 'CurlySMILES' and add the associated tags
        if (smiles[i] == '{' && (topology.size() != 0)) {
            auto prop_end = smiles.find_first_of('}', i);

            if (prop_end == std::string::npos) {
                throw format_error("SMI Reader: unmatched curly brace");
            }

            auto curly_prop = smiles.substr(i + 1, prop_end - i - 1);
            topology[topology.size() - 1].set("curly_property", std::string(curly_prop));
            i = prop_end;
            continue;
        }

        // We are not in a property list!
        // Therefore, if something is lowercase, then it must be
        // a single character element in the organic subset.
        if (is_ascii_lowercase(smiles[i])) {
            if (!is_aromatic_organic(smiles[i])) {
                throw format_error("SMI Reader: aromatic element not found: '{}'", smiles[i]);
            }
            auto name = smiles.substr(i, 1);
            this->add_atom(topology, name).set("is_aromatic", true);
            continue;
        }

        // We are not in a property list!
        // Therefore, aliphatic atoms can be written without brackets.
        // Br and Cl are tricky as, but the rest are single character
        if (is_ascii_uppercase(smiles[i])) {
            size_t element_length = 1;
            if (smiles[i] == 'C'
                && i + 1 < smiles.size()
                && smiles[i + 1] == 'l') {
                element_length = 2;
            }
            if (smiles[i] == 'B'
                && i + 1 < smiles.size()
                && smiles[i + 1] == 'r') {
                element_length = 2;
            }

            auto element_name = smiles.substr(i, element_length);
            if (!is_aliphatic_organic(element_name)) {
                throw format_error("SMI Reader: bare non-organic atom: '{}'", element_name);
            }

            this->add_atom(topology, element_name);
            i += element_length - 1;
            continue;
        }

        if (is_ascii_digit(smiles[i])) {
            auto ring_id = static_cast<size_t>(smiles[i] - '0');
            check_ring_(topology, ring_id);
            continue;
        }

        switch (smiles[i]) {
        case '.':
            if (!first_atom_) {
                first_atom_ = true;
                current_atom_++;
            }
            residues_.push_back(Residue("group " + std::to_string(groupid)));
            break;
        case '>': ++groupid; break;
        case '/': current_bond_order_ = Bond::UP; break;
        case '\\': current_bond_order_ = Bond::DOWN; break;
        case '~': current_bond_order_ = Bond::UNKNOWN; break;
        case '-':
            if (i + 1 < smiles.size() && smiles[i] == '>') {
                current_bond_order_ = Bond::DATIVE_R;
                ++i;
            } else {
                current_bond_order_ = Bond::SINGLE;
            }
            break;
        case '<':
            if (i + 1 < smiles.size() && smiles[i] == '-') {
                current_bond_order_ = Bond::DATIVE_L;
                ++i;
            } else {
                ++groupid; // Allow for reaction smiles, otherwise it's an unknown character
                first_atom_ = true;
            }
            break;
        case '=': current_bond_order_ = Bond::DOUBLE; break;
        case '#': current_bond_order_ = Bond::TRIPLE; break;
        case '$': current_bond_order_ = Bond::QUADRUPLE; break;
        case ':': current_bond_order_ = Bond::AROMATIC; break;
        case '(':
            branch_point_.push(previous_atom_);
            break;
        case ')':
            if (branch_point_.empty()) {
                throw format_error("SMI Reader: unmatched ')'");
            }
            previous_atom_ = branch_point_.top();
            branch_point_.pop();
            break;
        case '%':
            if (i + 2 >= smiles.size()) {
                throw format_error("SMI Reader: rings defined with '%' must be double digits");
            }
            check_ring_(topology, parse<size_t>(smiles.substr(i + 1, 2)));
            i += 2; // this line alone fixes most of issue 303 :)
            break;
        case '*':
            this->add_atom(topology, "*").set("wildcard", true);
            break;
        case ']':
        case '+':
        case '@':
            throw format_error("SMI Reader: symbol not allowed outside of property: '{}'", smiles[i]);
        default:
            throw format_error("SMI Reader: unknown symbol: '{}'", smiles[i]);
        }
    }

    for (auto residue: std::move(residues_)) {
        topology.add_residue(std::move(residue));
    }

    if (!branch_point_.empty()) {
        throw format_error("SMI Reader: {} unclosed '('(s)", branch_point_.size());
    }

    if (!rings_ids_.empty()) {
        throw format_error("SMI Reader: unclosed ring id '{}'", rings_ids_.begin()->first);
    }

    frame.resize(topology.size());
    frame.set_topology(topology);

    if (i < smiles.size()) {
        auto name = smiles.substr(i);
        frame.set("name", std::string(trim(name)));
    }
}

static void find_rings(
    const std::vector<std::vector<size_t>>& adj_list,
    std::unordered_map<size_t, size_t>& ring_atoms) {

    ring_atoms.clear();
    std::deque<bool> hit_atoms(adj_list.size(), false);
    std::set<Bond> ring_bonds;

    while (!all(hit_atoms)) {
        auto not_hit = std::find(hit_atoms.begin(), hit_atoms.end(), false);
        auto current_atom = static_cast<size_t>(std::distance(hit_atoms.begin(), not_hit));

        // mark this atom as processed
        hit_atoms[current_atom] = true;

        // it has no connections... therefore no rings!
        if (adj_list[current_atom].empty()) {
            continue;
        }

        // Search to find all rings for the component (all connected atoms)
        // for the 'current' atom.
        // Not guaranteed to fast, nor efficient, but just to work
        // Also not guaranteed to find the SSSR, but it will find all rings
        // needed for a given structure to be written as SMILES
        // adj_list The adjacency list graph
        // hit_atoms Atoms already encountered
        // ring_bonds Bonds already known to ring forming bonds
        // ring_atoms Map from an atom to the number of rings it forms

        // The atoms to process consist of :
        // current_atom Atom to be checked in this iteration
        // previous_atom Atom which was bound to the current atom

        std::stack<std::pair<size_t, size_t>> atoms_to_process;
        atoms_to_process.emplace(current_atom, current_atom);

        while (!atoms_to_process.empty()) {

            size_t previous_atom;

            std::tie(current_atom, previous_atom) = atoms_to_process.top();
            atoms_to_process.pop();

            const auto& current_atom_bonds = adj_list[current_atom];

            // We go through the neighbors backwards because we are using a
            // stack and want to go through them in the forward direction
            // when the stack is popped
            // (think of it like the negative of a negative is a positive)
            // We cannot use a queue as we want to process all neighbors of an
            // atom first, then move on to the next atom group (allowing for a
            // depth first like search instead of a breadth first search).
            for (auto neighbor_iter = current_atom_bonds.rbegin();
                neighbor_iter != current_atom_bonds.rend();
                ++neighbor_iter) {

                auto neighbor = *neighbor_iter;

                // prevent back tracking
                if (neighbor == previous_atom) {
                    continue;
                }

                // We've seen this neighbor before! Ring found
                // but only if we have NOT processed the current atom as this
                // prevents rings from being created from multiple directions.
                if (hit_atoms[neighbor] && !hit_atoms[current_atom]) {

                    // Don't process the same ring connection twice
                    if (ring_bonds.count(Bond(neighbor, current_atom)) != 0) {
                        continue;
                    }
                    ring_bonds.insert(Bond(neighbor, current_atom));

                    auto neigh_ring_iter = ring_atoms.find(neighbor);
                    if (neigh_ring_iter == ring_atoms.end()) {
                        neigh_ring_iter = ring_atoms.insert({ neighbor, 0 }).first;
                    }
                    neigh_ring_iter->second++;

                    continue; // No need to see this atom again
                }

                // Continue the search for new atoms
                if (!hit_atoms[neighbor]) {
                    // Add this atom to be processed. Replaces recursive section
                    atoms_to_process.push({ neighbor, current_atom });
                }
            }

            // We have processed this atom, mark it as done
            hit_atoms[current_atom] = true;
        }
    }
}

/// Sets int_part to `value` if `value` is actually an integer
static bool double_to_int(double value, int& int_part) {
    double integer_section;
    auto remainder = std::modf(value, &integer_section);
    if (remainder == 0.0) {
        int_part = static_cast<int>(std::floor(integer_section));
        return true;
    }
    return false;
}

static void write_atom_smiles(TextFile& file, const Atom& atom) {
    bool needs_brackets = false;
    auto type = atom.type();

    // The mass must be an integer is the only check we need as all atoms in the
    // periodic table have non-integer masses. Therefore, if the the mass is
    // an integer, then we know the user has set an isotope.
    int mass_int = 0;
    if (double_to_int(atom.mass(), mass_int)) {
        needs_brackets = true;
    }

    // If any of these values are set / not a default value
    auto chirality = atom.get<Property::STRING>("chirality").value_or("");
    if (!chirality.empty()) {
        needs_brackets = true;
    }

    int smi_class = -1;
    auto smi_class_prop = atom.get("smiles_class");
    if (smi_class_prop && smi_class_prop->kind() == Property::DOUBLE &&
        double_to_int(smi_class_prop->as_double(), smi_class) && smi_class >= 0) {
        needs_brackets = true;
    } else if (smi_class_prop) {
        warning("SMI Writer", "the 'smiles_class' property must be an integer >= 0");
    }

    int explicit_h = -1;
    auto explicit_h_prop = atom.get("hydrogen_count");
    if (explicit_h_prop && explicit_h_prop->kind() == Property::DOUBLE &&
        double_to_int(explicit_h_prop->as_double(), explicit_h) && explicit_h >= 0) {
        needs_brackets = true;
    } else if (explicit_h_prop) {
        warning("SMI Writer", "the 'hydrogen_count' property must be an integer >= 0");
    }

    // Charge is similar to mass. It must be an integer, otherwise it is ignored
    int charge_int = 0;
    if (double_to_int(atom.charge(), charge_int)) {
        // Use |= because we don't want to unset if charge_int is 0
        needs_brackets |= (charge_int != 0);
    }

    // Before the atom is printed, we must know if we are printing a property set.
    if (!is_aliphatic_organic(type)) {
        needs_brackets = true;
    }

    auto is_aromatic = atom.get<Property::BOOL>("is_aromatic").value_or(false);

    if (is_aromatic) {
        // We need to print brackets for aromatic Te, Se, As, and Si
        // Note, we allow aromatic Boron to be bare (no brackets) in following
        // ChemAxon's standard.
        if (type.size() > 1) {
            needs_brackets = true;
        }
        to_ascii_lowercase(type);
    }

    if (needs_brackets) {
        file.print("[");
    }

    // Mass must be first, before the element is printed
    if (mass_int != 0) {
        file.print("{}", mass_int);
    }

    file.print("{}", type);

    if (smi_class != -1) {
        file.print(":{}", smi_class);
    }

    bool is_good_tag = false;
    switch(chirality.size()) {
    case 0:
        is_good_tag = true;
        break;
    case 2:
        // If it's clockwise, then we just print the symbol @@ and be done with it
        if (chirality == "CW") {
            is_good_tag = true;
            file.print("@@");
        }
        break;
    case 3:
        if (chirality == "CCW") {
            is_good_tag = true;
            file.print("@");
        }
        break;
    case 7:
        if (chirality.find("CCW") == 0
            && is_chirality_tag(chirality.substr(4, 2))
            && is_ascii_digit(chirality[6])) {
            is_good_tag = true;
            file.print("@{}", chirality.substr(4));
        }
        break;
    case 8:
        if (chirality.find("CCW") == 0
            && is_chirality_tag(chirality.substr(4, 2))
            && is_ascii_digit(chirality[6])
            && is_ascii_digit(chirality[7])) {
            is_good_tag = true;
            file.print("@{}", chirality.substr(4));
        }
        break;
    default:
        break;
    }

    if (!is_good_tag) {
        warning("SMI Writer", "invalid chirality tag '{}'", chirality);
    }

    if (explicit_h == 1) {
        file.print("H");
    }

    if (explicit_h > 1) {
        file.print("H{}", explicit_h);
    }

    if (charge_int < 0) {
        file.print("-", charge_int);
    }

    if (charge_int > 0) {
        file.print("+", charge_int);
    }

    charge_int = std::abs(charge_int);
    if (charge_int != 1 && charge_int != 0) {
        file.print("{}", charge_int);
    }

    if (needs_brackets) {
        file.print("]");
    }
}

static void print_bond(TextFile& file, chemfiles::Bond::BondOrder bo) {
    switch(bo) {
    case Bond::SINGLE:
    case Bond::AMIDE:
        break;
    case Bond::DOUBLE: file.print("="); break;
    case Bond::TRIPLE: file.print("#"); break;
    case Bond::QUADRUPLE: file.print("$"); break;
    case Bond::AROMATIC: file.print(":"); break;
    case Bond::DATIVE_L: file.print("<-"); break;
    case Bond::DATIVE_R: file.print("->"); break;
    case Bond::UP: file.print("/"); break;
    case Bond::DOWN: file.print("\\"); break;
    case Bond::UNKNOWN:
        file.print("~");
        break;
    default:
        warning("SMI Writer", "unknown bond type");
        file.print("~");
        break;
    }
}

void SMIFormat::write_next(const Frame& frame) {
    if (frame.size() == 0) {
        file_.print("\n");
        return;
    }

    // Create an adjacent list as working with this is easier
    adj_list_.clear();
    adj_list_.resize(frame.size());
    for (const auto& bond : frame.topology().bonds()) {
        adj_list_[bond[0]].push_back(bond[1]);
        adj_list_[bond[1]].push_back(bond[0]);
    }

    find_rings(adj_list_, ring_atoms_);

    std::deque<bool> written(frame.size(), false);
    size_t branch_stack = 0;

    ring_stack_.clear();
    ring_count_ = 0;

    first_atom_= true;

    while (!all(written)) {
        if (!first_atom_) {
            file_.print(".");
        }
        auto not_hit = std::find(written.begin(), written.end(), false);
        auto start_atom = static_cast<size_t>(std::distance(written.begin(), not_hit));

        // We have found an atom that has not yet been printed! Now we must print it out
        // along with all of its connections (the entire component).

        // A structure to store a depth first like search through the component.
        std::stack<std::tuple<size_t, size_t, bool>> atoms_to_process;
        atoms_to_process.push(std::tuple<size_t,size_t,bool>( start_atom, start_atom, false ));

        while (!atoms_to_process.empty()) {
            bool needs_branch; // Do we need the '(' character?
            size_t previous_atom; // The atom that added the current atom to the stack
            size_t current_atom; // The current atom to process

            std::tie(previous_atom, current_atom, needs_branch) = atoms_to_process.top();
            atoms_to_process.pop();

            if (written[current_atom]) {
                continue;
            }

            auto& current_atom_bonds = adj_list_[current_atom];
            written[current_atom] = true;

            if (needs_branch) {
                file_.print("(");
                branch_stack++;
            }

            // This is only true the first time the loop is run
            if (current_atom != previous_atom) {
                print_bond(file_, frame.topology().bond_order(previous_atom, current_atom));
            }
            write_atom_smiles(file_, frame[current_atom]);

            // Prevent printing of additional '('
            size_t ring_start = 0;

            // We already know all rings in the structure, so the number of potential rings
            // for the current atom can be printed.
            auto any_rings = ring_atoms_.find(current_atom);
            if (any_rings != ring_atoms_.end()) {
                for (size_t i = 0; i < any_rings->second; ++i) {
                    ring_count_++;
                    ring_start++;
                    file_.print("{}", ring_count_ >= 10 ?
                                "%" + std::to_string(ring_count_) : std::to_string(ring_count_));
                    ring_stack_.insert({ current_atom, ring_count_ });
                }
            }

            // Avoid the printing of branch begin/end
            size_t ring_end = 0;

            // Find all ring connections first
            for (auto neighbor : current_atom_bonds) {
                // avoid 'trivial' rings
                if (neighbor == previous_atom) {
                    continue;
                }

                // We must have a ring to terminate
                if (written[neighbor]) {
                    auto ring = ring_stack_.find(neighbor);
                    if (ring != ring_stack_.end()) {
                        print_bond(file_,
                            frame.topology().bond_order(current_atom, neighbor)
                        );
                        file_.print("{}", ring->second >= 10 ?
                            "%" + std::to_string(ring->second) : std::to_string(ring->second));
                        ring_stack_.erase(ring);
                        ring_end++;
                    }
                }
            }

            // Handle branching
            size_t neighbors_printed = 0;
            for (auto neighbor_iter = current_atom_bonds.rbegin();
                neighbor_iter != current_atom_bonds.rend();
                ++neighbor_iter) {

                auto neighbor = *neighbor_iter;

                // prevent back tracking
                if (neighbor == previous_atom) {
                    continue;
                }

                // This got taken care of by printing a ring
                if (written[neighbor]) {
                    continue;
                }

                // To print a start bracket, we need to be branching (> 2 non-ring bonds)
                // and we don't want to branch the last neighbor printed
                auto needs_to_branch = neighbors_printed != 0 && neighbors_printed > ring_start;

                // Depth First Search like recursion
                atoms_to_process.push(std::tuple<size_t,size_t,bool>( current_atom, neighbor, needs_to_branch));

                // we printed a neighbor, if there's more than 1 neighbor, then we need to
                // branch for all but the last neighbor
                neighbors_printed++;
            }

            // End of branch
            if (current_atom_bonds.size() - ring_end == 1 && branch_stack != 0) {
                file_.print(")");
                branch_stack--;
            }
        }

        first_atom_ = false;
    }

    auto name = frame.get("name");
    if (name && name->kind() == Property::STRING) {
        file_.print("\t{}", name->as_string());
    }

    file_.print("\n");
}

optional<uint64_t> SMIFormat::forward() {
    auto position = file_.tellpos();

    auto line = file_.readline();
    while (trim(line).empty()) {
        if (file_.eof()) {
            return nullopt;
        }
        line = file_.readline();
    }

    return position;
}

bool all(const std::deque<bool>& vec) {
    for (auto i : vec) {
        if (!i) {
            return false;
        }
    }

    return true;
}
