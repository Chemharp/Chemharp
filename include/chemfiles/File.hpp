// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_FILES_HPP
#define CHEMFILES_FILES_HPP

#include <cstddef>
#include <cstdint>

#include <string>
#include <utility>
#include <vector>
#include <memory>
#include <string_view>

#include <fmt/format.h>

#include "chemfiles/exports.h"

namespace chemfiles {

/// Abstract base class for file representation.
class CHFL_EXPORT File {
public:
    /// Possible modes for opening a file
    enum Mode: char {
        /// Open in read-only mode
        READ = 'r',
        /// Open in read-write mode, and replace the file if it is already present
        WRITE = 'w',
        /// Open in read-write mode, and append new data at the end of the file
        APPEND = 'a',
    };

    /// Possible compression methods for opening a file
    enum Compression {
        /// Default method: plain text or binary formats
        DEFAULT,
        /// gzip compression
        GZIP,
        /// bzip2 compression
        BZIP2,
        /// lzma compression (.xz)
        LZMA,
    };

    virtual ~File() noexcept = default;
    File(File&&) = default;
    File& operator=(File&&) = default;
    File(const File&) = delete;
    File& operator=(const File&) = delete;

    /// Get the file path used to open this file.
    const std::string& path() const { return path_; }

    /// Get the mode used to open this file.
    Mode mode() const { return mode_; }

    /// Get the compression used to open this file.
    Compression compression() const { return compression_; }

protected:
    File(std::string path, Mode mode, Compression compression):
        path_(std::move(path)), mode_(mode), compression_(compression) {}

private:
    std::string path_;
    Mode mode_;
    Compression compression_;
};

/// Abstract base class for readers used by text files. This is specialized for
/// compressed files, and might get extended to network or memory mapped files.
///
/// Any failed operation should throw a `chemfiles::FileError`.
class TextFileImpl {
public:
    TextFileImpl(const std::string& path): path_(path) {}
    virtual ~TextFileImpl() = default;
    TextFileImpl(TextFileImpl&&) = delete;
    TextFileImpl(const TextFileImpl&) = delete;
    TextFileImpl& operator=(TextFileImpl&&) = delete;
    TextFileImpl& operator=(const TextFileImpl&) = delete;

    /// clear error and eof flags on the file
    virtual void clear() noexcept = 0;

    /// Set file indicator to `position` characters after the start of the file.
    ///
    /// @throws FileError in case of I/O error while seeking
    virtual void seek(uint64_t position) = 0;

    /// Fill up the `data` buffer reading at most `count` characters from the
    /// file.
    ///
    /// Returns the amount of characters read. A return value of 0 indicates
    /// that the end of file has been reached.
    ///
    /// @throws FileError in case of I/O error
    virtual size_t read(char* data, size_t count) = 0;

    /// Write `count` characters starting at `data` to the file.
    ///
    /// @throws FileError if it could not write all of the data to the file
    virtual void write(const char* data, size_t count) = 0;

protected:
    /// Get the string path used to open this file
    std::string_view path() const {
        return path_;
    }

private:
    std::string path_;
};

class MemoryBuffer;

/// Line-oriented text file reader and writer, using buffered read and fast
/// lines search.
///
/// This class reads text file line by line. It does so in a efficient way by
/// storing a chunk of the file in a memory buffer, and searching for new line
/// indicators (either `\n` or `\r\n`) in this buffer. It then returns
/// `string_view` inside this buffer, removing the need to allocate a new
/// `std::string` for each line.
///
/// Writing to the files is done without buffering or considering lines.
///
///
/// This class can read compressed data or in-memory data by using one of the
/// `TextFileImpl` interface implementation.
class CHFL_EXPORT TextFile final: public File {
public:
    /// Open the file at the given `path` with the requested `mode` and
    /// `compression` method, picking the best `TextFileImpl`.
    ///
    /// @throws FileError if the file does not exist in read mode, or if the
    ///                   file is invalid for the given compression method
    TextFile(std::string path, File::Mode mode, File::Compression compression);

    /// Use the given `MemoryBuffer` with the requested `mode` and `compression`
    /// method. A `MemoryFile` will be used as the `TextFileImpl`.
    ///
    /// @throws FileError if the file mode is append, or if trying to write to a
    /// compressed file.
    TextFile(std::shared_ptr<MemoryBuffer> memory, File::Mode mode, File::Compression compression);

    TextFile(TextFile&&) = default;
    TextFile& operator=(TextFile&&) = default;
    TextFile(const TextFile&) = delete;
    TextFile& operator=(const TextFile&) = delete;

    /// Returns the current position indicator, i.e. the number of characters
    /// from the beginning of the file.
    uint64_t tellpos() const;

    /// Set the position indicator to `position`.
    void seekpos(uint64_t position);

    /// Reset the position indicator to the beginning of the file, and clear
    /// end-of-file flag.
    void rewind();

    /// Check if end-of-file has been reached.
    bool eof() const {
        return eof_;
    }

    /// Clear end-of-file flags on the file.
    void clear();

    /// Read a single line from the file. The returned `string_view` points into
    /// an internal buffer, and can be invalidated after another call to
    /// `readline`. If storing the line is necessary, transform it to an owned
    /// string using `string_view::to_string()`.
    std::string_view readline();

    /// Read the full file into an owned string. This is a convenience method
    /// for format that need the full file read before parsing can start.
    std::string readall();

    /// Format and write some data to the file using the `fmt` library.
    ///
    /// This function has the same interface as `fmt::print(...)`, writting data
    /// to the file instead of stdout.
    template <typename Str, typename... Args>
    void print(const Str& format, const Args&... args) {
        this->vprint(format, fmt::make_format_args(args...));
    }

private:
    /// Fill the buffer, calling `refill` and setting all needed internal values
    void fill_buffer(size_t start);

    /// Actually format and print data to the file
    void vprint(fmt::string_view format, fmt::format_args args);

    /// Check if the buffer was initialized and contains data from the
    /// underlying  `TextFileImpl`.
    bool buffer_initialized() const;

    /// Pointer to the actual file implementation
    std::unique_ptr<TextFileImpl> file_;
    /// Buffer storing characters read from the `TextFileImpl`. If
    /// `got_impl_eof_` is true, this contains the remaining characters from the
    /// file and is then padded with null characters ('\0').
    std::vector<char> buffer_;
    /// Start of the current line;
    const char* line_start_;
    /// End of the buffer
    const char* end_;
    /// Current position in the file, this is the number of characters in the
    /// file up to the start of the buffer (`buffer_[0]` is the `position_`'th
    /// character in the file).
    uint64_t position_ = 0;
    /// Did we reach the end of the underlying `TextFileImpl`? Since we are
    /// buffering data, this does not necessarily correspond to `this->eof()`.
    bool got_impl_eof_ = false;
    /// Did we actually reached the end of file while reading a line?
    bool eof_ = false;
};

} // namespace chemfiles

#endif
