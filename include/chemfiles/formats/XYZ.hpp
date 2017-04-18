// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_FORMAT_XYZ_HPP
#define CHEMFILES_FORMAT_XYZ_HPP

#include <string>

#include "chemfiles/Format.hpp"
#include "chemfiles/FormatFactory.hpp"

namespace chemfiles {

/// [XYZ] file format reader and writer.
///
/// [XYZ]: http://openbabel.org/wiki/XYZ
class XYZFormat final: public Format {
public:
    XYZFormat(const std::string& path, File::Mode mode);

    void read_step(size_t step, Frame& frame) override;
    void read(Frame& frame) override;
    void write(const Frame& frame) override;
    std::string description() const override;
    size_t nsteps() override;

    // Register the xyz format with the ".xyz" extension and the "XYZ" description.
    FORMAT_NAME(XYZ)
    FORMAT_EXTENSION(.xyz)
private:
    /// Text file where we read from
    std::unique_ptr<TextFile> file_;
    /// Storing the positions of all the steps in the file, so that we can
    /// just `seekg` them instead of reading the whole step.
    std::vector<std::streampos> steps_positions_;
};

} // namespace chemfiles

#endif
