// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <string>
#include <string_view>

#include "chemfiles/FormatFactory.hpp"
#include "chemfiles/FormatMetadata.hpp"
#include "chemfiles/File.hpp"

#include "chemfiles/misc.hpp"
#include "chemfiles/error_fmt.hpp"
#include "chemfiles/external/optional.hpp"

using namespace chemfiles;

/// try to distinguish CIF and mmCIF files, since they share the same
/// `.cif` extension
static optional<std::string> distinguish_cif_variants(const std::string& path, const std::string& compression);

static bool contains(std::string_view haystack, std::string_view needle) {
    return haystack.find(needle) != std::string_view::npos;
}

std::string chemfiles::guess_format(std::string path, char mode) {
    std::string extension;
    std::string compression;

    auto dot1 = path.rfind('.');
    if (dot1 != std::string::npos) {
        extension = path.substr(dot1);
        bool new_extension = false;
        // check file extension for compressed file extension
        if (extension == ".gz") {
            new_extension = true;
            compression = "GZ";
        } else if (extension == ".bz2") {
            new_extension = true;
            compression = "BZ2";
        } else if (extension == ".xz") {
            new_extension = true;
            compression = "XZ";
        }

        if (new_extension) {
            extension = "";
            auto dot2 = path.substr(0, dot1).rfind('.');
            if (dot2 != std::string::npos) {
                extension = path.substr(0, dot1).substr(dot2);
            }
        }
    }

    if (extension.empty()) {
        throw file_error(
            "file at '{}' does not have an extension, provide a format name to read it",
            path
        );
    }

    if (extension == ".cif" && (mode == 'r' || mode == 'a')) {
        extension = distinguish_cif_variants(path, compression).value_or(extension);
    }

    auto registered_format = FormatFactory::get().by_extension(extension);
    auto format = std::string(registered_format.metadata.name);

    if (!compression.empty()) {
        format += " / " + compression;
    }

    return format;
}


static optional<std::string> distinguish_cif_variants(
    const std::string& path,
    const std::string& compression
) {
    auto c = File::DEFAULT;
    if (compression == "GZ") {
        c = File::GZIP;
    } else if (compression == "BZ2") {
        c = File::BZIP2;
    } else if (compression == "XZ") {
        c = File::LZMA;
    }

    try {
        auto file = TextFile(path, File::READ, c);
        while (!file.eof()) {
            auto line = file.readline();
            // check a few mmCIF/CIF specific tags that are more likely to be
            // close to the top of the file
            if (contains(line, "_audit_conform.dict_name")
                || contains(line, "_cell.length_a")
                || contains(line, "_atom_site.type_symbol")
            ) {
                return std::string(".mmcif");
            }

            if (contains(line, "_symmetry_equiv_pos_as_xyz")
                || contains(line, "_cell_length_a")
                || contains(line, "_atom_site_type_symbol")
            ) {
                return std::string(".cif");
            }
        }

        // if we could not find any of the above strings in the file, it is very
        // likely the file is invalid. As below, the user will get a proper
        // error when trying to open the file
        return nullopt;
    } catch (const FileError&) {
        // in case of error while reading, just use the file extension for now,
        // the user will get a proper error when trying to open the file
        return nullopt;
    }
}
