// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_FORMAT_CIF_HPP
#define CHEMFILES_FORMAT_CIF_HPP

#include "chemfiles/config.h"  // IWYU pragma: keep
#ifndef CHFL_DISABLE_GEMMI

#include <cstddef>
#include <string>
#include <vector>
#include <memory>
#include <utility>

#include "gemmi/small.hpp"

#include "chemfiles/File.hpp"
#include "chemfiles/Format.hpp"

namespace chemfiles {
class Frame;
class MemoryBuffer;
class FormatMetadata;

/// CIF (Crystallographic Information Framework) files reader and writer.
///
/// The reader code is based on the [gemmi](https://project-gemmi.github.io/)
/// project.
class CIFFormat final: public Format {
public:
    CIFFormat(std::string path, File::Mode mode, File::Compression compression) :
        file_(std::move(path), mode, compression), current_step_(0) {
        init_();
    }

    CIFFormat(std::shared_ptr<MemoryBuffer> memory, File::Mode mode, File::Compression compression) :
        file_(std::move(memory), mode, compression), current_step_(0) {
        init_();
    }

    void read_step(size_t step, Frame& frame) override;
    void read(Frame& frame) override;
    void write(const Frame& frame) override;
    size_t nsteps() override;
private:
    /// Initialize important variables
    void init_();
    /// Underlying file representation
    TextFile file_;
    /// Store all structures, reading the whole file during init_()
    std::vector<gemmi::SmallStructure> structures_;
    /// When reading frame by frame, or writing, remember where we are
    size_t current_step_;
};

template<> const FormatMetadata& format_metadata<CIFFormat>();

} // namespace chemfiles

#endif // CHFL_DISABLE_GEMMI

#endif
