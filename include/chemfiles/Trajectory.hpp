// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_TRAJECTORY_HPP
#define CHEMFILES_TRAJECTORY_HPP

#include <cstddef>
#include <memory>
#include <string>

#include "chemfiles/exports.h"
#include "chemfiles/Frame.hpp"
#include "chemfiles/UnitCell.hpp"
#include "chemfiles/external/span.hpp"  // IWYU pragma: keep
#include "chemfiles/external/optional.hpp"

namespace chemfiles {
class Format;
class Topology;
class MemoryBuffer;

/// A `Trajectory` is a chemistry file on the hard drive. It is the entry point
/// of the chemfiles library.
class CHFL_EXPORT Trajectory final {
public:
    /// Open a file, automatically gessing the file format from the extension.
    ///
    /// The `format` parameter should be a string formatted as `"<format>"` or
    /// `"<format>/<compression>"`. `<format>` should be the format name (see
    /// the corresponding [documentation section][formats] for the names) or an
    /// empty string. `<compression>` should be `GZ` for gzip files, `BZ2` for
    /// bzip2 files, or `XZ` for lzma/.xz files. If `<compression>` is present,
    /// it will determine which compression method is used to read/write the
    /// file.
    ///
    /// For example, `format = "XYZ"` will force usage of XYZ format regardless
    /// of the file extension; `format = "XYZ / GZ"` will additionally use gzip
    /// compression.
    ///
    /// If the `format` is an empty string, the file extension will be used to
    /// guess the format. If the file path ends with `.gz`, `.xz`, or `.bz2`;
    /// the file will be treated as a compressed file and the next extension is
    /// used to guess the format. For example `Trajectory("file.xyz.gz")` will
    /// open the file for reading using the XYZ format and the gzip compression
    /// method.
    ///
    /// @example{trajectory/trajectory.cpp}
    ///
    /// @param path The file path. In `w` or `a` modes, the file is
    ///             created if it does not exist yet. In `r` mode, an
    ///             exception is thrown is the file does not exist yet.
    /// @param mode Opening mode for the file. Default mode is `r` for read.
    ///             Other supported modes depends on the underlying format and
    ///             are `w` for write, and `a` for append. `w` mode discard any
    ///             previously existing file.
    /// @param format Specific format to use. Needed when there is no way to
    ///               guess the format from the extension of the file, or when
    ///               this guess would be wrong.
    ///
    /// @throws FileError for all errors concerning the physical file: can not
    ///                   open it, can not read/write it, *etc.*
    /// @throws FormatError if the file is not valid for the used format.
    ///
    /// [formats]: http://chemfiles.org/chemfiles/latest/formats.html#list-of-supported-formats
    explicit Trajectory(std::string path, char mode = 'r', const std::string& format = "");

    /// Read a memory buffer as though it were a formatted file
    ///
    /// The `format` parameter should be follow the same rules as in the main
    /// `Trajectory` constructor.
    ///
    /// @example{trajectory/memory_reader.cpp}
    ///
    /// @param data The start of the memory buffer used to store the file.
    ///             It does not need to be *null* terminated.
    /// @param size The size of the memory buffer.
    /// @param format Specific format to use.
    ///
    /// @throws FileError If the compression given in `format` is not supported
    /// @throws FormatError if the data in the buffer is not valid for the used
    ///                     format, or the format does not support reading from
    ///                     a memory buffer
    static Trajectory memory_reader(const char* data, size_t size, const std::string& format);

    /// Write to a memory buffer as though it were a formatted file
    ///
    /// The `format` parameter should be follow the same rules as in the main
    /// `Trajectory` constructor, except that compression specification are not
    /// supported.
    ///
    /// To retreive the memory written to by the returned `Trajectory` object,
    /// make a call to the `memory_buffer` function.
    ///
    /// @example{trajectory/memory_writer.cpp}
    ///
    /// @param format Specific format to use.
    ///
    /// @throws FileError If any compression is given in `format`
    /// @throws FormatError if the format does not support writing to a memory buffer
    static Trajectory memory_writer(const std::string& format);

    ~Trajectory();

    Trajectory(Trajectory&& other) noexcept;
    Trajectory& operator=(Trajectory&& other) noexcept;

    Trajectory(const Trajectory& other) = delete;
    Trajectory& operator=(const Trajectory& other) = delete;

    /// Read the next frame in the trajectory.
    ///
    /// The trajectory must have been opened in read mode, and the
    /// underlying format must support reading.
    ///
    /// This function throws a `FileError` if there are no more frames to read
    /// in the trajectory.
    ///
    /// @example{trajectory/read.cpp}
    ///
    /// @throws FileError for all errors concerning the physical file: can not
    ///                   open it, can not read/write it, *etc.*
    /// @throws FormatError if the file is not valid for the used format, or if
    ///                     the format does not support reading.
    Frame read();

    /// Read a single frame at specified `step` from the trajectory.
    ///
    /// The trajectory must have been opened in read mode, and the
    /// underlying format must support reading.
    ///
    /// This function throws a `FileError` if the step is bigger than the
    /// number of steps in the trajectory.
    ///
    /// @example{trajectory/read_step.cpp}
    ///
    /// @param step step to read from the trajectory
    ///
    /// @throws FileError for all errors concerning the physical file: can not
    ///                   open it, can not read/write it, *etc.*
    /// @throws FormatError if the file is not valid for the used format, or if
    ///                     the format does not support reading.
    Frame read_step(size_t step);

    /// Write a single frame to the trajectory.
    ///
    /// The trajectory must have been opened in write or append mode, and the
    /// underlying format must support reading.
    ///
    /// @example{trajectory/write.cpp}
    ///
    /// @param frame frame to write to this trajectory
    ///
    /// @throws FileError for all errors concerning the physical file: can not
    ///                   open it, can not read/write it, *etc.*
    /// @throws FormatError if the format does not support writing.
    void write(const Frame& frame);

    /// Use the given `topology` instead of any pre-existing `Topology` when
    /// reading or writing.
    ///
    /// This replace any topology in the file being read, or in the `Frame`
    /// being written.
    ///
    /// This is mainly usefull when a format does not define topological
    /// information, as it can be the case with some molecular dynamic formats.
    ///
    /// @example{trajectory/set_topology.cpp}
    ///
    /// @param topology the topology to use with this frame
    ///
    /// @throws Error if the topology does not contain the right number of
    ///               atoms at any step.
    void set_topology(const Topology& topology);

    /// Use the `Topology` from the first `Frame` of the `Trajectory` at
    /// `filename` instead any pre-existing `Topology` when reading or writing.
    ///
    /// This replace any topology in the file being read, or in the `Frame`
    /// being written.
    ///
    /// This is mainly usefull when a format does not define topological
    /// information, as it can be the case with some molecular dynamic formats.
    ///
    /// @example{trajectory/set_topology.cpp}
    ///
    /// @param filename trajectory file path.
    /// @param format Specific format to use. Needed when there is no way to
    ///               guess the format from the extension of the file, or when
    ///               this guess would be wrong.
    ///
    /// @throws FileError for all errors concerning the physical file: can not
    ///                   open it, can not read it, *etc.*
    /// @throws FormatError if the file is not valid for the used format.
    /// @throws Error if the topology does not contain the right number of
    ///               atoms at any step.
    void set_topology(const std::string& filename, const std::string& format = "");

    /// Use the given `cell` instead of any pre-existing `UnitCell` when
    /// reading or writing.
    ///
    /// This replace any unit cell in the file being read, or in the `Frame`
    /// being written.
    ///
    /// This is mainly usefull when a format does not define unit cell
    /// information.
    ///
    /// @param cell the unit cell to use with this frame
    ///
    /// @example{trajectory/set_cell.cpp}
    void set_cell(const UnitCell& cell);

    /// Get the number of steps (the number of frames) in this trajectory.
    ///
    /// @example{trajectory/nsteps.cpp}
    size_t nsteps() const;

    /// Check if all the frames in this trajectory have been read, *i.e.* if
    /// the last read frame is the last frame of the trajectory.
    ///
    /// @example{trajectory/done.cpp}
    bool done() const;

    /// Close a trajectory, and synchronize all buffered content with the drive.
    ///
    /// Calling any function on a closed trajectory will throw a `FileError`.
    ///
    /// @example{trajectory/close.cpp}
    void close();

    /// Get the path used to open the trajectory
    ///
    /// @example{trajectory/path.cpp}
    const std::string& path() const {
        return path_;
    }

    /// Get the memory buffer used for writing if the trajectory was created
    /// with `Trajectory::memory_writer`.
    ///
    /// If the trajectory was not created for writing to memory, this will
    /// return `nullopt`.
    ///
    /// @example{trajectory/memory_buffer.cpp}
    optional<span<const char>> memory_buffer() const;

private:
    Trajectory(char mode, std::unique_ptr<Format> format, std::shared_ptr<MemoryBuffer> buffer);

    /// Perform a few checks before reading a frame
    void pre_read(size_t step);
    /// Set the frame topology and/or cell after reading it
    void post_read(Frame& frame);
    /// Check that the trajectory is still open, and throw a `FileError` is it
    /// has been closed.
    void check_opened() const;

    /// Path of the associated file
    std::string path_;
    /// Opening mode of the associated file
    char mode_ = '\0';
    /// Current step
    size_t step_ = 0;
    /// Number of steps in the file, if available
    size_t nsteps_ = 0;
    /// Format used to read the associated file. It will be `nullptr` is the
    /// trajectory is closed
    std::unique_ptr<Format> format_;
    /// Topology to use for reading/writing files when no topological data is
    /// present
    optional<Topology> custom_topology_;
    /// UnitCell to use for reading/writing files when no unit cell information
    /// is present
    optional<UnitCell> custom_cell_;
    /// The internal memory buffer, shared with the MemoryFile implementation
    std::shared_ptr<MemoryBuffer> buffer_;
};

} // namespace chemfiles

#endif
