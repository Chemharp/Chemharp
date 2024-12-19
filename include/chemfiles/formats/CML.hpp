// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_FORMAT_CML_HPP
#define CHEMFILES_FORMAT_CML_HPP

#include <cstddef>
#include <string>
#include <memory>
#include <unordered_map>
#include <utility>

#include <pugixml.hpp>

#include "chemfiles/File.hpp"
#include "chemfiles/Format.hpp"

namespace chemfiles {
class Frame;
class MemoryBuffer;
class FormatMetadata;

/// CML file format reader.
class CMLFormat final: public Format {
public:
    CMLFormat(std::string path, File::Mode mode, File::Compression compression):
        file_(std::move(path), mode, compression) {
        init_();
    }

    CMLFormat(std::shared_ptr<MemoryBuffer> memory, File::Mode mode, File::Compression compression):
        file_(std::move(memory), mode, compression) {
        init_();
    }

    ~CMLFormat() override;

    void read_step(size_t step, Frame& frame) override;
    void read(Frame& frame) override;
    void write(const Frame& frame) override;
    size_t nsteps() override;
private:
    /// Initialize the document and root objects
    void init_();

    /// Read the atoms from `atoms` into `frame`
    void read_atoms(Frame& frame, const pugi::xml_node& atoms);
    /// Read the bonds from `bonds` into `frame`
    void read_bonds(Frame& frame, const pugi::xml_node& bonds);

    /// Text file where we read from. It needs to stay valid if we write to the file
    TextFile file_;

    /// XML document must be kept alive with the object
    pugi::xml_document document_;

    /// First actual node. Could be 'cml' (multi-frame) or 'molecule'(single-frame).
    /// Others (like single atom) are not supported
    pugi::xml_node root_;

    /// If multi-frame, store the current location
    pugi::xml_named_node_iterator current_;

    /// Number of frames added to the file
    size_t num_added_ = 0;

    /// Atomic reference to atomic id in the current frame
    std::unordered_map<std::string, size_t> ref_to_id_;
};

template<> const FormatMetadata& format_metadata<CMLFormat>();

} // namespace chemfiles

#endif
