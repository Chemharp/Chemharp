// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_FORMAT_HPP
#define CHEMFILES_FORMAT_HPP

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>
#include <memory>

#include "chemfiles/exports.h"

#include "chemfiles/File.hpp"
#include "chemfiles/Error.hpp"
#include "chemfiles/external/optional.hpp"

namespace chemfiles {
class Frame;
class MemoryBuffer;
class FormatMetadata;

/// Get the metadata associated with `Format`.
///
/// The metadata should be a reference to static memory.
///
/// In order to implement a new format, one should specialise this function
/// with the corresponding format:
///
/// ```cpp
/// class MyFormat: public Format {
///     // ...
/// };
///
/// namespace chemfiles {
///     template<> const FormatMetadata& format_metadata<MyFormat>() {
///         static FormatMetadata metadata;
///         metadata.name = "MyFormat";
///         metadata.extension = ".mtf";
///         return metadata;
///     }
/// }
/// ```
template<class Format>
const FormatMetadata& format_metadata() {
    throw FormatError("format_metadata is not implemented for this format");
}

/// The `Format` class defines the interface to implement in order to add a new
/// format to chemfiles. For text-based formats, it might be simpler to
/// implement the `TextFormat` interface instead.
///
/// It is possible to implement only one of `Format::read`; `Format::read_step`
/// or `Format::write`. In that case, only the corresponding operations will be
/// available from the corresponding `chemfiles::Trajectory`.
class CHFL_EXPORT Format {
public:
    Format() = default;
    virtual ~Format() = default;

    // Delete all move and copy constructors.
    // Format should only be used behind an std::unique_ptr
    Format(const Format&) = delete;
    Format& operator=(const Format&) = delete;
    Format(Format&&) = delete;
    Format& operator=(Format&&) = delete;

    /// Read a specific `step` from the trajectory file.
    ///
    /// @throw FormatError if the file does not follow the format
    /// @throw FileError if their is an OS error while reading the file
    ///
    /// @param step The step to read
    /// @param frame The frame to fill
    virtual void read_step(size_t step, Frame& frame);

    /// Read the next step from the trajectory file.
    ///
    /// @throw FormatError if the file does not follow the format
    /// @throw FileError if their is an OS error while reading the file
    ///
    /// @param frame The frame to fill
    virtual void read(Frame& frame);

    /// Write a frame to the trajectory file.
    ///
    /// @throw FormatError if the file does not follow the format
    /// @throw FileError if their is an OS error while reading the file
    ///
    /// @param frame The frame to be written
    virtual void write(const Frame& frame);

    /// Get the number of frames in the associated file. This function can be
    /// expensive to call since it may needs to scan the whole file.
    ///
    /// @return The number of frames
    virtual size_t nsteps() = 0;
};

/// The `TextFormat` class defines a common, simpler interface for text based
/// formats.
///
/// It is possible to implement only one of `TextFormat::read_next` or
/// `TextFormat::write_next`. In that case, only the corresponding operations
/// will be available from the corresponding `chemfiles::Trajectory`.
class CHFL_EXPORT TextFormat: public Format {
public:
    TextFormat(std::string path, File::Mode mode, File::Compression compression);
    TextFormat(std::shared_ptr<MemoryBuffer> memory, File::Mode mode, File::Compression compression);
    virtual ~TextFormat() override = default;

    void read_step(size_t step, Frame& frame) override;
    void read(Frame& frame) override;
    void write(const Frame& frame) override;
    size_t nsteps() override;

    /// Fast-forward the file for one step, returning a valid position if the
    /// file does contain one more step or `nullopt` if it does not.
    virtual optional<uint64_t> forward() = 0;

    virtual void read_next(Frame& frame);
    virtual void write_next(const Frame& frame);

protected:
    /// Text file used to read/write data
    TextFile file_;

private:
    /// Scan the whole file to get all the steps positions
    void scan_all();

    /// The next step to read
    size_t step_ = 0;

    /// Storing the positions of all the steps in the file, so that we can
    /// just `seekpos` them instead of reading the whole step.
    std::vector<uint64_t> steps_positions_;

    /// Did we found the end of file while scanning or reading?
    bool eof_found_ = false;
};

} // namespace chemfiles

#endif
