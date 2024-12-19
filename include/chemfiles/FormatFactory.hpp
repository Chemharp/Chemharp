// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_FORMAT_FACTORY_HPP
#define CHEMFILES_FORMAT_FACTORY_HPP

#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <functional>
#include <type_traits>

#include "chemfiles/exports.h"
#include "chemfiles/mutex.hpp"
#include "chemfiles/File.hpp"
#include "chemfiles/Format.hpp"

namespace chemfiles {
class MemoryBuffer;
class FormatMetadata;

using format_creator_t = std::function<std::unique_ptr<Format>(std::string path, File::Mode mode, File::Compression compression)>;
using memory_stream_t = std::function<std::unique_ptr<Format>(std::shared_ptr<MemoryBuffer> memory, File::Mode mode, File::Compression compression)>;

struct RegisteredFormat {
    const FormatMetadata& metadata;
    format_creator_t creator;
    memory_stream_t memory_stream_creator;
};

template <typename T>
using SupportsMemoryIO = std::is_constructible<T, std::shared_ptr<MemoryBuffer>, File::Mode, File::Compression>;

/// This class allow to register Format with names and file extensions
class CHFL_EXPORT FormatFactory final {
private:
    FormatFactory();

public:
    /// Get the instance of the `FormatFactory`
    static FormatFactory& get();

    /// Get a `RegisteredFormat` from a format name.
    ///
    /// @param name the format name
    /// @throws FormatError if the format can not be found
    const RegisteredFormat& by_name(const std::string& name);

    /// Get a `RegisteredFormat` from a format extention.
    ///
    /// @param extension the format extention
    /// @throws FormatError if the format can not be found
    const RegisteredFormat& by_extension(const std::string& extension);

    /// Register a given `Format` in the factory if it supports memory IO
    ///
    /// The format informations are taken from the specialization of the
    /// template function `chemfiles::format_metadata` for this format. The
    /// second template argument is used to determine if the `Format` supports
    /// memory IO.
    template<class Format, std::enable_if_t<SupportsMemoryIO<Format>::value, int> = 0>
    void add_format() {
        const auto& metadata = format_metadata<Format>();
        metadata.validate();
        register_format(metadata,
            [](const std::string& path, File::Mode mode, File::Compression compression) {
                return std::make_unique<Format>(path, mode, compression);
            },
            [](std::shared_ptr<MemoryBuffer> memory, File::Mode mode, File::Compression compression) {
                return std::make_unique<Format>(std::move(memory), mode, compression);
            }
        );
    }

    /// Register a given `Format` in the factory if it does not support memory IO
    ///
    /// The format informations are taken from the specialization of the
    /// template function `chemfiles::format_metadata` for this format. The
    /// second template argument is used to determine if the `Format` supports
    /// memory IO.
    template<class Format, std::enable_if_t<!SupportsMemoryIO<Format>::value, int> = 0>
    void add_format() {
        const auto& metadata = format_metadata<Format>();
        metadata.validate();
        register_format(metadata,
            [](const std::string& path, File::Mode mode, File::Compression compression) {
                return std::make_unique<Format>(path, mode, compression);
            }
        );
    }

    /// Get the metadata for all registered formats
    std::vector<std::reference_wrapper<const FormatMetadata>> formats();

private:
    void register_format(const FormatMetadata& metadata, format_creator_t creator, memory_stream_t memory_stream);
    void register_format(const FormatMetadata& metadata, format_creator_t creator);

    /// Trajectory map associating format descriptions and creators
    mutex<std::vector<RegisteredFormat>> formats_;
};

} // namespace chemfiles

#endif
