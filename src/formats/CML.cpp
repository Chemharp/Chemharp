// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cmath>
#include <cstddef>
#include <array>
#include <string>
#include <utility>
#include <vector>
#include <iterator>
#include <string_view>
#include <unordered_map>

#include <pugixml.hpp>

#include "chemfiles/types.hpp"
#include "chemfiles/unreachable.hpp"
#include "chemfiles/utils.hpp"
#include "chemfiles/parse.hpp"
#include "chemfiles/warnings.hpp"
#include "chemfiles/error_fmt.hpp"
#include "chemfiles/external/optional.hpp"

#include "chemfiles/File.hpp"
#include "chemfiles/Atom.hpp"
#include "chemfiles/Frame.hpp"
#include "chemfiles/Property.hpp"
#include "chemfiles/Topology.hpp"
#include "chemfiles/UnitCell.hpp"
#include "chemfiles/Connectivity.hpp"
#include "chemfiles/FormatMetadata.hpp"
#include "chemfiles/periodic_table.hpp"

#include "chemfiles/formats/CML.hpp"

using namespace chemfiles;

template<> const FormatMetadata& chemfiles::format_metadata<CMLFormat>() {
    static FormatMetadata metadata;
    metadata.name = "CML";
    metadata.extension = ".cml";
    metadata.description = "Chemical Markup Language";
    metadata.reference = "http://www.xml-cml.org";

    metadata.read = true;
    metadata.write = true;
    metadata.memory = true;

    metadata.positions = true;
    metadata.velocities = true;
    metadata.unit_cell = true;
    metadata.atoms = true;
    metadata.bonds = true;
    metadata.residues = false;
    return metadata;
}

class xml_writer final: public pugi::xml_writer {
public:
    xml_writer(TextFile& file): file_(file) {}

    void write(const void* data, size_t size) override {
        file_.print("{}", std::string_view(static_cast<const char*>(data), size));
    }

private:
    TextFile& file_;
};

void CMLFormat::init_() {
    if (file_.mode() == File::WRITE) {
        root_ = document_.append_child("cml");
        root_.append_attribute("xmlns") = "http://www.xml-cml.org/schema";
        root_.append_attribute("xmlns:cml") = "http://www.xml-cml.org/dict/cml";
        root_.append_attribute("xmlns:units") = "http://www.xml-cml.org/units/units";
        root_.append_attribute("xmlns:convention") = "http://www.xml-cml.org/convention";
        root_.append_attribute("convention") = "convention:molecular";
        root_.append_attribute("xmlns:iupac") = "http://www.iupac.org";
        return;
    }

    // may need to check the file for a CML node in the future for a CML node,
    // and if it exists append to that node and rewrite the file, appending new
    // frames to the CML node. For now, let's just add molecule nodes to the end
    // of the file.
    if (file_.mode() == File::APPEND) {
        root_ = document_;
        return;
    }

    auto content = file_.readall();

    auto result = document_.load_string(content.c_str());
    if (!result) {
        throw format_error("parsing error: '{}'", result.description());
    }

    root_ = document_.child("cml");
    if (root_ != nullptr) {
        auto molecules = root_.children("molecule");
        current_ = molecules.begin();
        if (current_ == molecules.end()) {
            throw format_error("cml node has no valid children");
        }
        return;
    }

    // Not technically standard to have multiple molecule nodes as root nodes,
    // but some tools do this and we should support it
    if (document_.child("molecule") != nullptr) {
        auto molecules = document_.children("molecule");
        current_ = molecules.begin();
        root_ = document_;
        return;
    }

    throw format_error("no supported starting node found");
}

CMLFormat::~CMLFormat() {
    if (num_added_ != 0 && file_.mode() == File::WRITE) {
        auto writter = xml_writer(file_);
        document_.save(writter, "  ");
    } else if (file_.mode() == File::APPEND) {
        // Don't bother to check if anything is added, the document will be blank
        // regardless
        auto writter = xml_writer(file_);
        document_.save(writter, "  ", pugi::format_no_declaration | pugi::format_default);
    }
}

size_t CMLFormat::nsteps() {
    auto children = root_.children("molecule");
    return static_cast<size_t>(std::distance(children.begin(), children.end()));
}

// Read a property into either a frame or an atom with supplied title
template<typename T>
static void read_property_(T& container, const pugi::xml_node& node, const std::string& title) {
    std::string data_type_str = "xsd:string";
    auto data_type = node.attribute("dataType");
    if (!data_type) {
        warning("CML reader", "{} has no data type, assuming string", title);
    } else {
        data_type_str = data_type.as_string();
    }

    auto value = node.text();
    if (data_type_str == "xsd:string") {
        container.set(title, value.as_string());
    } else if (data_type_str == "xsd:boolean") {
        container.set(title, value.as_bool());
    } else if (data_type_str == "xsd:double" || data_type_str == "xsd:integer") {
        container.set(title, value.as_double());
    } else {
        warning("CML reader", "unknown data type: {}, assuming string", data_type_str);
        container.set(title, value.as_string());
    }
}

// Read a property into either a frame or an atom
template<typename T>
static void read_property_(T& container, const pugi::xml_node& node) {
    auto title = node.attribute("title");
    if (!title) {
        warning("CML reader", "skipping untitled property");
        return;
    }

    read_property_(container, node, title.as_string());
}

static UnitCell read_cell(const pugi::xml_node& crystal) {
    Vector3D lengths;
    Vector3D angles = {90, 90, 90};
    for (const auto& scalar: crystal.children("scalar")) {
        auto title = scalar.attribute("title");
        if (title != nullptr) {
            auto name = std::string(title.value());
            if (name == "a") {
                lengths[0] = scalar.text().as_double();
            } else if (name == "b") {
                lengths[1] = scalar.text().as_double();
            } else if (name == "c") {
                lengths[2] = scalar.text().as_double();
            } else if (name == "alpha") {
                angles[0] = scalar.text().as_double();
            } else if (name == "beta") {
                angles[1] = scalar.text().as_double();
            } else if (name == "gamma") {
                angles[2] = scalar.text().as_double();
            } else {
                warning("CML reader", "unknown crystal scalar: {}", name);
            }
        }
    }
    return UnitCell(lengths, angles);
}

void CMLFormat::read_atoms(Frame& frame, const pugi::xml_node& atoms) {
    ref_to_id_.clear();

    // TODO: Read string vector versions instead of just the nodes
    for (const auto& atom: atoms.children("atom")) {
        double x2 = 0;
        double y2 = 0;
        double x3 = 0;
        double y3 = 0;
        double z3 = 0;
        double xf = 0;
        double yf = 0;
        double zf = 0;
        std::string id;
        std::string element;
        std::string atom_title;
        double formal_charge = 0.0;
        double isotope = 0.0;
        double hydrogen_count = 0.0;
        for (const auto& attribute: atom.attributes()) {
            auto name = std::string(attribute.name());
            if (name == "id") {
                id = attribute.as_string();
            } else if (name == "elementType") {
                element = attribute.as_string();
            } else if (name == "x2") {
                x2 = attribute.as_double();
            } else if (name == "y2") {
                y2 = attribute.as_double();
            } else if (name == "x3") {
                x3 = attribute.as_double();
            } else if (name == "y3") {
                y3 = attribute.as_double();
            } else if (name == "z3") {
                z3 = attribute.as_double();
            } else if (name == "xFract") {
                xf = attribute.as_double();
            } else if (name == "yFract") {
                yf = attribute.as_double();
            } else if (name == "zFract") {
                zf = attribute.as_double();
            } else if (name == "title") { // special attributes
                atom_title = attribute.as_string();
            } else if (name == "formalCharge") {
                formal_charge = attribute.as_double();
            } else if (name == "isotopeNumber") {
                isotope = attribute.as_double();
            } else if (name == "hydrogenCount") {
                hydrogen_count = attribute.as_double();
            } else {
                warning("CML reader", "unknown atom attribute: {}", name);
            }
        }

        auto new_atom = Atom(atom_title, element);

        Vector3D position;
        if (frame.cell().shape() != UnitCell::INFINITE && (xf != 0.0 || yf != 0.0 || zf != 0.0)) {
            position = frame.cell().matrix() * Vector3D(xf, yf, zf);
        } else if (x3 == 0.0 && y3 == 0.0 && z3 == 0.0) {
            position = Vector3D(x2, y2, 0);
        } else {
            position = Vector3D(x3, y3, z3);
        }

        // Set special attributes
        if (formal_charge != 0.0) {
            new_atom.set_charge(formal_charge);
        }

        if (isotope != 0.0) {
            new_atom.set_mass(isotope);
        }

        if (hydrogen_count != 0.0) {
            new_atom.set("hydrogen_count", hydrogen_count);
        }

        for (const auto& scalar: atom.children("scalar")) {
            read_property_(new_atom, scalar);
        }

        // Vectors are stored separate where the values of the vector are
        // stored as doubles separated by whitespace.
        Vector3D velocity;
        for (const auto& vector3: atom.children("vector3")) {
            const auto& title_attribute = vector3.attribute("title");
            if (!title_attribute) {
                warning("CML reader", "vector3 is not titled");
                continue;
            }
            std::string title = title_attribute.as_string();

            auto vect_strings = split(vector3.text().as_string(), ' ');
            if (vect_strings.size() != 3) {
                warning("CML reader", "{} vector3 does not have 3 values", title);
                continue;
            }

            Vector3D vector;
            try {
                vector[0] = parse<double>(vect_strings[0]);
                vector[1] = parse<double>(vect_strings[1]);
                vector[2] = parse<double>(vect_strings[2]);
            } catch (const Error&) {
                warning("CML reader", "{} contains elements which are not numeric", title);
                continue;
            }

            if (title == "velocity") {
                velocity = vector;
                frame.add_velocities();
            } else {
                new_atom.set(title, vector);
            }
        }

        frame.add_atom(std::move(new_atom), position, velocity);
        ref_to_id_[id] = frame.size() - 1;
    }
}

void CMLFormat::read_bonds(Frame& frame, const pugi::xml_node& bonds) {
    for (const auto& bond: bonds.children("bond")) {
        auto atomref = bond.attribute("atomRefs2");
        auto order = bond.attribute("order");

        if (!atomref) {
            warning("CML reader", "bond does not contain an atomref attribute");
            continue;
        }

        auto ids = split(atomref.as_string(), ' ');
        if (ids.size() != 2) {
            warning("CML reader", "bondArray contains a bond of size {} instead of 2", ids.size());
            continue;
        }

        auto id1 = ref_to_id_.find(std::string(ids[0]));
        auto id2 = ref_to_id_.find(std::string(ids[1]));
        if (id1 == ref_to_id_.end() || id2 == ref_to_id_.end()) {
            warning("CML reader", "invalid atomic references in bond: {} -- {}", ids[0], ids[1]);
            continue;
        }

        auto bond_order = Bond::UNKNOWN;
        if (order != nullptr) {
            std::string order_str = order.as_string();
            if (!order_str.empty()) {
                switch (order_str[0]) {
                case '1':
                case 'S':
                case 's':
                    bond_order = Bond::SINGLE;
                    break;
                case '2':
                case 'D':
                case 'd':
                    bond_order = Bond::DOUBLE;
                    break;
                case '3':
                case 'T':
                case 't':
                    bond_order = Bond::TRIPLE;
                    break;
                case 'A':
                case 'a':
                    bond_order = Bond::AROMATIC;
                    break;
                default:
                    warning("CML reader", "unknown bond order: '{}'", order_str[0]);
                    break;
                }
            }
        }

        frame.add_bond(id1->second, id2->second, bond_order);
    }
}

void CMLFormat::read(Frame& frame) {
    const auto& current = *current_;
    for (const auto& attribute: current.attributes()) {
        auto name = std::string(attribute.name());
        if (name == "id") { // do nothing, required junk for validation
        } else if (name == "xmlns") {
        } else if (name == "xmlns:cml") {
        } else if (name == "xmlns:units") {
        } else if (name == "xmlns:xsd") {
        } else if (name == "xmlns:iupac") {
        } else if (name == "title") {
            frame.set("title", attribute.as_string());
        } else {
            warning("CML reader", "unknown molecule attribute: {}", name);
        }
    }

    const auto& frame_name = current.child("name");
    if (frame_name != nullptr) {
        frame.set("name", frame_name.text().as_string());
    }

    const auto& crystal = current.child("crystal");
    if (crystal != nullptr) {
        frame.set_cell(read_cell(crystal));
    }

    // Frame level properties
    const auto& properties = current.child("propertyList");
    if (properties != nullptr) {
        for (const auto& property: properties.children("property")) {
            const auto& title = property.attribute("title");
            if (!title) {
                warning("CML reader", "skipping untitled property");
                continue;
            }

            const auto* title_str = title.as_string();
            const auto& scalar = property.child("scalar");
            if (!scalar) {
                warning("CML reader", "{} has no scalar associated with it", title_str);
                continue;
            }

            read_property_(frame, scalar, title_str);
        }
    }

    auto atoms = current.child("atomArray");
    if (atoms != nullptr) {
        read_atoms(frame, atoms);
    } else {
        warning("CML reader", "missing atomArray node");
    }

    auto bonds = current.child("bondArray");
    if (bonds != nullptr) {
        read_bonds(frame, bonds);
    }

    current_++;
}

void CMLFormat::read_step(size_t step, Frame& frame) {
    current_ = root_.children("molecule").begin();
    std::advance(current_, static_cast<ptrdiff_t>(step));
    read(frame);
}

static void write_property_(const Property& p, pugi::xml_node& node) {
    Vector3D v;
    switch (p.kind()) {
    case Property::DOUBLE:
        node.append_attribute("dataType") = "xsd:double";
        node.text() = p.as_double();
        break;
    case Property::BOOL:
        node.append_attribute("dataType") = "xsd:boolean";
        node.text() = p.as_bool();
        break;
    case Property::STRING:
        node.append_attribute("dataType") = "xsd:string";
        node.text() = p.as_string().c_str();
        break;
    case Property::VECTOR3D: // Supported only in schema3 for properties
        node.set_name("vector3");
        v = p.as_vector3d();
        node.text() = (std::to_string(v[0]) + " " +
                       std::to_string(v[1]) + " " +
                       std::to_string(v[2])).c_str();
        break;
    default:
        unreachable();
    }
}

// Used for isotopes and formal charges
static bool is_double_integer(double val) {
    return std::abs(val - std::trunc(val)) < 1e-3;
}

void CMLFormat::write(const Frame& frame) {
    auto mol = root_.append_child("molecule");

    if (file_.mode() == File::WRITE) {
        ++num_added_;
        mol.append_attribute("id") = ("m" + std::to_string(num_added_)).c_str();
    }

    // The title is supplied as a attribute
    auto title_prop = frame.get("title");
    if (title_prop && title_prop->kind() == Property::STRING) {
        mol.append_attribute("title") = title_prop->as_string().c_str();
    }

    // The name is supplied as a child node
    auto name_prop = frame.get("name");
    if (name_prop && name_prop->kind() == Property::STRING) {
        mol.append_child("name").text() = name_prop->as_string().c_str();
    }

    const auto& cell = frame.cell();
    if (cell.shape() != UnitCell::INFINITE) {
        auto lengths = cell.lengths();
        auto angles = cell.angles();

        auto crystal_node = mol.append_child("crystal");
        auto scalar = crystal_node.append_child("scalar");
        scalar.append_attribute("units") = "units:angstrom";
        scalar.append_attribute("title") = "a";
        scalar.text() = lengths[0];

        scalar = crystal_node.append_child("scalar");
        scalar.append_attribute("units") = "units:angstrom";
        scalar.append_attribute("title") = "b";
        scalar.text() = lengths[1];

        scalar = crystal_node.append_child("scalar");
        scalar.append_attribute("units") = "units:angstrom";
        scalar.append_attribute("title") = "c";
        scalar.text() = lengths[2];

        scalar = crystal_node.append_child("scalar");
        scalar.append_attribute("units") = "units:degree";
        scalar.append_attribute("title") = "alpha";
        scalar.text() = angles[0];

        scalar = crystal_node.append_child("scalar");
        scalar.append_attribute("units") = "units:degree";
        scalar.append_attribute("title") = "beta";
        scalar.text() = angles[1];

        scalar = crystal_node.append_child("scalar");
        scalar.append_attribute("units") = "units:degree";
        scalar.append_attribute("title") = "gamma";
        scalar.text() = angles[2];
    }

    // Loop through and add properties to the propertyList node as required
    const auto& properties = frame.properties();
    if (properties.size() != 0) {
        auto prop_list = mol.append_child("propertyList");
        for (const auto& prop: properties) {

            // This properties are special, don't write them as scalars
            if (prop.first == "name" || prop.first == "title") {
                continue;
            }

            auto prop_node = prop_list.append_child("property");
            prop_node.append_attribute("title") = prop.first.c_str();

            auto scalar_node = prop_node.append_child("scalar");
            write_property_(prop.second, scalar_node);
        }
    }

    auto velocities = frame.velocities();
    auto atom_array = mol.append_child("atomArray");
    for (size_t i = 0; i < frame.size(); ++i) {
        const auto& atom = frame[i];
        auto atom_node = atom_array.append_child("atom");

        // Add a new id for the atom which is 1-based
        atom_node.append_attribute("id") = ("a" + std::to_string(i + 1)).c_str();
        atom_node.append_attribute("elementType") = atom.type().c_str();

        // Write special attributes here instead of as scalars
        if (!atom.name().empty()) {
            atom_node.append_attribute("title") = atom.name().c_str();
        }

        // Make sure charge is set and is an integer
        if (atom.charge() != 0.0 && is_double_integer(atom.charge())) {
            atom_node.append_attribute("formalCharge") = static_cast<int>(atom.charge());
        }

        // See if the mass of the atom is set to something non-default and is an integer
        auto atom_info = find_in_periodic_table(atom.type());
        if (atom_info && atom_info->mass &&
            atom.mass() != *(atom_info->mass) &&
            is_double_integer(atom.mass())) {
            atom_node.append_attribute("isotopeNumber") = static_cast<size_t>(atom.mass());
        }

        auto hydrogen_count = atom.get("hydrogen_count");
        if (hydrogen_count && hydrogen_count->kind() == Property::DOUBLE) {
            atom_node.append_attribute("hydrogenCount") =
                static_cast<size_t>(hydrogen_count->as_double());
        }

        // Now write the position attributes
        const auto& position = frame.positions()[i];
        atom_node.append_attribute("x3") = position[0];
        atom_node.append_attribute("y3") = position[1];
        atom_node.append_attribute("z3") = position[2];

        if (velocities) { // write velocities if they are included
            auto node = atom_node.append_child("vector3");
            node.append_attribute("title") = "velocity";

            auto v = velocities->at(i);
            node.text() = (std::to_string(v[0]) + " " +
                           std::to_string(v[1]) + " " +
                           std::to_string(v[2])).c_str();

        }

        const auto& atom_properties = atom.properties();
        if (!atom_properties) {
            continue;
        }
        for (const auto& prop : *atom_properties) {
            // special properties which are already written as attributes
            // charge and isotope are not stored as properties, so no need to check
            if (prop.first == "hydrogen_count" || prop.first == "title") {
                continue;
            }
            auto scalar_node = atom_node.append_child("scalar");
            scalar_node.append_attribute("title") = prop.first.c_str();
            write_property_(prop.second, scalar_node);
        }
    }

    const auto& bonds = frame.topology().bonds();
    if (bonds.size() == 0) { // don't bother if there's no bonds
        return;
    }

    auto bond_array = mol.append_child("bondArray");
    const auto& bond_orders = frame.topology().bond_orders();

    for (size_t i = 0; i < bonds.size(); ++i) {
        std::string refs2 = "a" + std::to_string(bonds[i][0] + 1) +
                           " a" + std::to_string(bonds[i][1] + 1);
        auto bond_node = bond_array.append_child("bond");
        bond_node.append_attribute("atomRefs2") = refs2.c_str();

        switch (bond_orders[i]) {
            case Bond::SINGLE:
                bond_node.append_attribute("order") = "1";
                break;
            case Bond::DOUBLE:
                bond_node.append_attribute("order") = "2";
                break;
            case Bond::TRIPLE:
                bond_node.append_attribute("order") = "3";
                break;
            case Bond::AROMATIC:
                bond_node.append_attribute("order") = "a";
                break;
            default:
                break;
        }
    }
}
