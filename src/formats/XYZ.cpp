// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdlib>

#include <map>
#include <set>
#include <array>
#include <string>
#include <utility>
#include <vector>
#include <string_view>
#include <unordered_map>

#include <fmt/core.h>
#include <fmt/format.h>

#include "chemfiles/types.hpp"
#include "chemfiles/parse.hpp"
#include "chemfiles/utils.hpp"
#include "chemfiles/warnings.hpp"
#include "chemfiles/error_fmt.hpp"
#include "chemfiles/unreachable.hpp"
#include "chemfiles/external/optional.hpp"

#include "chemfiles/File.hpp"
#include "chemfiles/Atom.hpp"
#include "chemfiles/Frame.hpp"
#include "chemfiles/Property.hpp"
#include "chemfiles/UnitCell.hpp"
#include "chemfiles/FormatMetadata.hpp"

#include "chemfiles/formats/XYZ.hpp"

using namespace chemfiles;

/// Parse the Lattice="..." in extended XYZ
static UnitCell parse_cell(std::string_view lattice);

struct extended_property {
    extended_property(std::string name_, Property::Kind type_): name(std::move(name_)), type(type_) {}
    std::string name;
    Property::Kind type;
};

/// Type for a list of additional atomic properties in extended XYZ
using properties_list_t = std::vector<extended_property>;

/// Read the extended XYZ properties from the line, set frame properties
/// directly and return the list of atomic properties to read
static properties_list_t read_extended_comment_line(std::string_view line, Frame& frame);

/// Read the properties in the list form the line and set them on the atom
static void read_atomic_properties(const properties_list_t& properties, std::string_view line, Atom& atom);

/// Get the list of atoms properties defined for all atoms in the frame
static properties_list_t get_atom_properties(const Frame& frame);

/// Generate the extended XYZ comment line for the given frame
static std::string write_extended_comment_line(const Frame& frame, const properties_list_t& properties);

template<> const FormatMetadata& chemfiles::format_metadata<XYZFormat>() {
    static FormatMetadata metadata;
    metadata.name = "XYZ";
    metadata.extension = ".xyz";
    metadata.description = "XYZ text format";
    metadata.reference = "https://openbabel.org/wiki/XYZ";

    metadata.read = true;
    metadata.write = true;
    metadata.memory = true;

    metadata.positions = true;
    metadata.velocities = false;
    metadata.unit_cell = true;
    metadata.atoms = true;
    metadata.bonds = false;
    metadata.residues = false;
    return metadata;
}

void XYZFormat::read_next(Frame& frame) {
    auto n_atoms = parse<size_t>(file_.readline());

    auto properties = read_extended_comment_line(file_.readline(), frame);

    frame.reserve(n_atoms);
    for (size_t i=0; i<n_atoms; i++) {
        auto line = file_.readline();
        double x = 0;
        double y = 0;
        double z = 0;
        std::string name;
        auto count = scan(line, name, x, y, z);
        auto atom = Atom(std::move(name));
        read_atomic_properties(properties, line.substr(count), atom);
        frame.add_atom(std::move(atom), Vector3D(x, y, z));
    }
}

void XYZFormat::write_next(const Frame& frame) {
    const auto& positions = frame.positions();
    auto properties = get_atom_properties(frame);

    file_.print("{}\n", frame.size());
    file_.print("{}\n", write_extended_comment_line(frame, properties));

    for (size_t i = 0; i < frame.size(); i++) {
        const auto& atom = frame[i];

        auto name = atom.name();
        if (name.empty()) {
            name = "X";
        }

        file_.print("{} {:g} {:g} {:g}",
            name, positions[i][0], positions[i][1], positions[i][2]
        );

        for (const auto& property: properties) {
            const auto& value = atom.get(property.name).value();

            if (property.type == Property::STRING) {
                file_.print(" {}", value.as_string());
            } else if (property.type == Property::BOOL) {
                if (value.as_bool()) {
                    file_.print(" T");
                } else {
                    file_.print(" F");
                }
            } else if (property.type == Property::DOUBLE) {
                file_.print(" {:g}", value.as_double());
            } else if (property.type == Property::VECTOR3D) {
                const auto& vector = value.as_vector3d();
                file_.print(" {:g} {:g} {:g}", vector[0], vector[1], vector[2]);
            }
        }

        file_.print("\n");
    }
}

optional<uint64_t> XYZFormat::forward() {
    auto position = file_.tellpos();

    auto line = file_.readline();
    if (trim(line).empty() || file_.eof()) {
        // We just read an empty line, give up here
        return nullopt;
    }

    size_t n_atoms = 0;
    try {
        n_atoms = parse<size_t>(line);
    } catch (const Error&) {
        throw format_error(
            "could not read the number of atoms for XYZ format: the line is '{}'",
            line
        );
    }

    for (size_t i=0; i<n_atoms + 1; i++) {
        if (file_.eof()) {
            throw format_error(
                "XYZ format: not enough lines at step {} (expected {}, got {})",
                current_forward_step_, n_atoms + 2, i + 1
            );
        }

        file_.readline();
    }

    current_forward_step_++;

    return position;
}

/*****************************************************************************/
/* End of the main XYZ format implementation. The code below implements the  */
/* extended XYZ format/convention. This format is 100% backward compatible   */
/* with normal XYZ. It adds data in two places: after the atomic positions   */
/* and in the comment line. The comment line contains key/values pairs       */
/* describing frame-level properties, as well as two special keys: Lattice   */
/* stores the unit cell matrix, and 'Properties' describe the name and type  */
/* of any additional atom-level data.                                        */
/*****************************************************************************/

// Should we surround `s` with quotation marks when using it as property name?
static bool should_be_quoted(std::string_view s) {
    for (auto c: s) {
        // TODO: ASE also allow [] {} and () to function as quotes. This should
        // be updated when a specification is agreed on.
        if (is_ascii_whitespace(c) || c == '=' || c == '\'' || c == '"') {
            return true;
        }
    }
    return false;
}

static bool contains_single_quote(std::string_view s) {
    for (auto c: s) {
        if (c == '\'') {
            return true;
        }
    }
    return false;
}

static bool contains_double_quote(std::string_view s) {
    for (auto c: s) {
        if (c == '"') {
            return true;
        }
    }
    return false;
}

std::string write_extended_comment_line(const Frame& frame, const properties_list_t& properties) {
    std::string result = "Properties=species:S:1:pos:R:3";

    for (const auto& property: properties) {
        char type;
        int count;
        switch (property.type) {
        case Property::STRING:
            type = 'S';
            count = 1;
            break;
        case Property::BOOL:
            type = 'L';
            count = 1;
            break;
        case Property::DOUBLE:
            type = 'R';
            count = 1;
            break;
        case Property::VECTOR3D:
            type = 'R';
            count = 3;
            break;
        default:
            unreachable();
        }

        result += fmt::format(":{}:{}:{}", property.name, type, count);
    }

    // support for Lattice
    if (frame.cell().shape() != UnitCell::INFINITE) {
        auto matrix = frame.cell().matrix();
        // set small elements to 0
        for (size_t i=0; i<3; i++) {
            for (size_t j=0; j<3; j++) {
                if (std::abs(matrix[i][j]) < 1e-12) {
                    matrix[i][j] = 0;
                }
            }
        }
        result += fmt::format(" Lattice=\"{:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g}\"",
            matrix[0][0], matrix[1][0], matrix[2][0],
            matrix[0][1], matrix[1][1], matrix[2][1],
            matrix[0][2], matrix[1][2], matrix[2][2]
        );
    }

    // sort properties to have reproducible output
    auto sorted_properties = std::map<std::string_view, Property>{
        frame.properties().begin(), frame.properties().end()
    };
    // support for generic frame properties
    for (const auto& it: sorted_properties) {
        if (should_be_quoted(it.first)) {
            // quote the string
            if (!contains_double_quote(it.first)) {
                result += fmt::format(" \"{}\"=", it.first);
            } else if (!contains_single_quote(it.first)) {
                result += fmt::format(" '{}'=", it.first);
            } else {
                 warning(
                    "Extended XYZ",
                    "frame property '{}' contains both single and double quote, it will not be saved",
                    it.first
                );
                 continue;
            }
        } else {
            result += fmt::format(" {}=", it.first);
        }

        switch (it.second.kind()) {
        case Property::STRING:
            result += fmt::format("\"{}\"", it.second.as_string());
            break;
        case Property::BOOL:
            result += fmt::format("{}", it.second.as_bool() ? 'T' : 'F');
            break;
        case Property::DOUBLE:
            result += fmt::format("{:g}", it.second.as_double());
            break;
        case Property::VECTOR3D:
            {
                auto v = it.second.as_vector3d();
                result += fmt::format("\"{:g} {:g} {:g}\"", v[0], v[1], v[2]);
            }
            break;
        default:
            unreachable();
        }
    }
    return result;
}

static bool is_valid_property_name(const std::string& name) {
    if (name.empty()) {
        return false;
    }

    if (!is_ascii_letter(name[0])) {
        return false;
    }

    for (auto c: name) {
        if (!is_ascii_alphanumeric(c) && c != '_') {
            return false;
        }
    }

    return true;
}

properties_list_t get_atom_properties(const Frame& frame) {
    if (frame.size() == 0) {
        return {};
    }

    auto all_properties = std::map<std::string, Property::Kind>();
    auto partially_defined_already_warned = std::set<std::string>();

    const auto& first_atom = frame[0];
    if (first_atom.properties()) {
        for (const auto& property: *first_atom.properties()) {
            if (!is_valid_property_name(property.first)) {
                warning(
                    "Extended XYZ", "'{}' is not a valid property name for extended "
                    "XYZ, is will not be saved",
                    property.first, 0
                );
                partially_defined_already_warned.insert(property.first);
                continue;
            }

            if (property.second.kind() == Property::STRING) {
                if (should_be_quoted(property.second.as_string())) {
                    warning(
                        "Extended XYZ", "string value for property '{}' on atom {} "
                        "can not be be saved as an atomic property",
                        property.first, 0
                    );
                    continue;
                }
            }

            all_properties.emplace(property.first, property.second.kind());
        }
    }

    for (size_t i=1; i<frame.size(); i++) {
        const auto& atom = frame[i];

        auto to_remove = std::vector<std::string>();
        for (const auto& property: all_properties) {
            auto current_property = atom.get(property.first);
            if (!current_property) {
                // this property was present for all atoms until now, but
                // not on the current one
                warning(
                    "Extended XYZ",
                    "property '{}' is only defined for a subset of atoms, it will not be saved",
                    property.first
                );

                to_remove.push_back(property.first);
                continue;
            }

            if (current_property.value().kind() != property.second) {
                warning(
                    "Extended XYZ",
                    "property '{}' is defined with different types on different atoms, if will not be saved",
                    property.first
                );
                partially_defined_already_warned.insert(property.first);
                to_remove.push_back(property.first);
                continue;
            }
        }

        const auto& atom_properties = atom.properties();
        if (atom_properties && (*atom_properties).size() > all_properties.size()) {
            // warn for properties defined on this atom but not on others
            for (const auto& property: *atom.properties()) {
                if (all_properties.count(property.first) == 0) {
                    if (partially_defined_already_warned.count(property.first) == 0) {
                        warning(
                            "Extended XYZ",
                            "property '{}' is only defined for a subset of atoms, it will not be saved",
                            property.first
                        );
                        partially_defined_already_warned.insert(property.first);
                    }
                }
            }
        }

        for (const auto& name: to_remove) {
            all_properties.erase(name);
        }
    }

    auto results = properties_list_t();
    results.reserve(all_properties.size());
    for (auto property: std::move(all_properties)) {
        results.emplace_back(property.first, property.second);
    }
    return results;
}

/// mapping name=>properties values. The use of string_view is safe since the
/// original string_view (as returned by file.readline()) will stay alive longer
/// than this map.
using extended_xyz_properties_map_t = std::unordered_map<std::string_view, Property>;

/// A simple parser for the extended XYZ comment line. This parser bails out as
/// soon as possible if the line does not seems to follow extended XYZ
/// convention.
///
/// This line should contains key/values pairs separated by spaces, while the
/// key and the value are separated by an = sign. If the value contains spaces,
/// it is then quoted, both single quotes and double quotes are accepted here.
/// Example: Properties=species:S:1:pos:R:3:bad:I:1 Lattice="3 0 0 0 4 0 0 0 5" name='simple file'
struct extended_xyz_parser {
public:
    extended_xyz_parser(std::string_view line): current_(line.data()), end_(line.data() + line.size()) {}
    ~extended_xyz_parser() = default;
    extended_xyz_parser(extended_xyz_parser&&) = delete;
    extended_xyz_parser(const extended_xyz_parser&) = delete;
    extended_xyz_parser& operator=(extended_xyz_parser&&) = delete;
    extended_xyz_parser& operator=(const extended_xyz_parser&) = delete;

    extended_xyz_properties_map_t parse() {
        extended_xyz_properties_map_t properties;
        while (!done()) {
            skip_whitespace();
            auto name = next_substring();
            if (!name.empty() && current() == '=') {
                advance();
            } else {
                // properties without associated value are set to True by ASE
                properties.emplace(name, true);
                continue;
            }

            auto value = next_substring();
            if (!is_ascii_whitespace(current()) && !done()) {
                // missing a space after the value, bail out here
                warning("Extended XYZ", "expected whitespace after the value for {}, got {}", name, current());
                break;
            }

            properties.emplace(name, std::string(value));
        }
        return properties;
    }

private:
    void skip_whitespace() {
        while (!done()) {
            if (is_ascii_whitespace(current())) {
                advance();
            } else {
                return;
            }
        }
    }

    /// Get a string, either quoted and possibly containing spaces; or unquoted.
    /// Unquoted strings ends on a space or '=' sign.
    std::string_view next_substring() {
        bool check_for_quote = false;
        char quote = '\0';
        if (current() == '"' || current() == '\'') {
            check_for_quote = true;
            quote = current();
            advance();
        }

        const auto* start = current_;
        size_t size = 0;
        while (!done()) {
            if (check_for_quote && current() == quote) {
                // end of quoted value
                check_for_quote = false;
                advance();
                break;
            } else if (!check_for_quote && (is_ascii_whitespace(current()) || current() == '=')) {
                // end of non-quoted value
                break;
            } else {
                advance();
                size++;
            }
        }

        auto value = std::string_view(start, size);
        if (check_for_quote) {
            // we found an initial quote, but no final one
            assert(done());
            warning("Extended XYZ", "missing final quote after {}", value);
        }
        return value;
    }

    // are we done with the input?
    bool done() const {
        return current_ == end_;
    }

    // Advance the pointer and return the new character
    void advance() {
        if (!done()) {
            current_++;
        }
    }

    char current() {
        if (done()) {
            return '\0';
        } else {
            return *current_;
        }
    }

    const char* current_;
    const char* end_;
};


static UnitCell parse_cell(std::string_view lattice) {
    auto matrix = Matrix3D::zero();
    scan(lattice,
        matrix[0][0], matrix[1][0], matrix[2][0],
        matrix[0][1], matrix[1][1], matrix[2][1],
        matrix[0][2], matrix[1][2], matrix[2][2]
    );
    return UnitCell(matrix);
}

/// Parse extended XYZ properties format string
///
/// Format is "[NAME:TYPE:NCOLS]...", e.g. "species:S:1:pos:R:3:CS:R:4".
/// NAME is the name of the property.
/// TYPE is one of R, I, S, L for real, integer, string and logical.
/// NCOLS is number of columns for that property.
///
/// The standard XYZ properties ("species:S:1:pos:R:3".) are checked, but not
/// returned in the property list.
///
/// When converting atomic properties as described by the 'Properties'
/// value, the following rules are used:
///
/// - Properties of type L:1 are maped to boolean values
/// - Properties of type L:N are maped to N separate boolean properties,
///   named `$name_$i`
///
/// - Properties of type R:1 and I:1 are maped to double values
/// - Properties of type R:3 and I:3 are maped to Vector3D values
/// - Properties of type R:N and I:N are maped to to N separate double
///   properties, named `$name_$i`
///
/// - Properties of type S:1 are maped to string values
/// - Properties of type S:N are maped to N separate string properties,
///   named `$name_$i`
static properties_list_t parse_property_list(std::string_view input) {
    if (input.substr(0, 19) != "species:S:1:pos:R:3") {
        warning("Extended XYZ", "ignoring non-standard Properties='{}', should start with 'species:S:1:pos:R:3'", input);
        return {};
    }

    // saving the initial input for error messages
    auto list = input;
    if (list.size() <= 20) {
        return {};
    }
    list.remove_prefix(20);
    properties_list_t properties;

    auto splitted = split(list, ':');
    if (splitted.size() % 3 != 0) {
        warning("Extended XYZ", "ignoring invalid Properties='{}'", input);
        return {};
    }

    auto count = splitted.size() / 3;
    for (size_t i=0; i<count; i++) {
        auto name = std::string(splitted[3 * i]);

        auto type_str = splitted[3 * i + 1];
        Property::Kind type;
        if (type_str == "R" || type_str == "I") {
            type = Property::DOUBLE;
        } else if (type_str == "S") {
            type = Property::STRING;
        } else if (type_str == "L") {
            type = Property::BOOL;
        } else {
            warning("Extended XYZ", "invalid type name for {} in Properties='{}'", name, input);
            continue;
        }

        size_t repeat = 0;
        try {
            repeat = parse<size_t>(splitted[3 * i + 2]);
        } catch (const Error&) {
            warning("Extended XYZ", "invalid type repeat for {} in Properties='{}'", name, input);
            continue;
        }

        if (repeat == 3 && type == Property::DOUBLE) {
            type = Property::VECTOR3D;
            properties.emplace_back(extended_property{name, type});
            continue;
        }

        if (repeat == 0) {
            warning("Extended XYZ", "invalid type repeat for {} in Properties='{}'", name, input);
            continue;
        } else if (repeat == 1) {
            properties.emplace_back(extended_property{name, type});
        } else {
            for (size_t j=0; j<repeat; j++) {
                properties.emplace_back(extended_property{
                    fmt::format("{}_{}", name, j), type
                });
            }
        }
    }

    return properties;
}

properties_list_t read_extended_comment_line(std::string_view line, Frame& frame) {
    // only try to parse as extended XYZ if `Properties` is defined as expected
    if (line.find("species:S:1:pos:R:3") == std::string::npos) {
        return {};
    }
    auto properties = extended_xyz_parser(line).parse();

    for (const auto& it: properties) {
        auto name = it.first;
        if (name == "Lattice" || name == "Properties") {
            continue;
        }
        frame.set(std::string(name), it.second);
    }

    if (properties.count("Lattice") == 1) {
        frame.set_cell(parse_cell(properties.at("Lattice").as_string()));
    }

    if (properties.count("Properties") == 1) {
        return parse_property_list(properties.at("Properties").as_string());
    } else {
        return {};
    }
}

// This function throws `chemfiles::Error` if the values are not convertible to
// the expected type. If the files contains a valid `Properties=...`
// description,throwing errors if the rest of the files does not follow the
// description is fair game.
void read_atomic_properties(const properties_list_t& properties, std::string_view line, Atom& atom) {
    for (const auto& property: properties) {
        if (property.type == Property::STRING) {
            std::string value;
            auto count = scan(line, value);
            line.remove_prefix(count);
            atom.set(property.name, std::move(value));
        } else if (property.type == Property::BOOL) {
            std::string value;
            auto count = scan(line, value);
            line.remove_prefix(count);
            to_ascii_lowercase(value);
            if (value == "t" || value == "true") {
                atom.set(property.name, true);
            } else if (value == "f" || value == "false") {
                atom.set(property.name, false);
            } else {
                throw error("invalid value for boolean '{}'", value);
            }
        } else if (property.type == Property::DOUBLE) {
            double value;
            auto count = scan(line, value);
            line.remove_prefix(count);
            atom.set(property.name, value);
        }  else if (property.type == Property::VECTOR3D) {
            Vector3D value;
            auto count = scan(line, value[0], value[1], value[2]);
            line.remove_prefix(count);
            atom.set(property.name, value);
        } else {
            unreachable();
        }
    }
}
