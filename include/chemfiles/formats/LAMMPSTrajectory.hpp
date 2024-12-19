// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_FORMAT_LAMMPS_TRAJECTORY_HPP
#define CHEMFILES_FORMAT_LAMMPS_TRAJECTORY_HPP

#include <cstddef>
#include <cstdint>

#include <array>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>

#include "chemfiles/File.hpp"
#include "chemfiles/Format.hpp"

#include "chemfiles/external/optional.hpp"

namespace chemfiles {
class Frame;
class MemoryBuffer;
class FormatMetadata;

/// LAMMPS Atom file format reader and writer.
class LAMMPSTrajectoryFormat final : public TextFormat {
  public:
    LAMMPSTrajectoryFormat(std::string path, File::Mode mode, File::Compression compression)
        : TextFormat(std::move(path), mode, compression) {}

    LAMMPSTrajectoryFormat(std::shared_ptr<MemoryBuffer> memory, File::Mode mode,
                     File::Compression compression)
        : TextFormat(std::move(memory), mode, compression) {}

    void read_next(Frame& frame) override;
    void write_next(const Frame& frame) override;
    optional<uint64_t> forward() override;

  private:
    std::array<double, 3> read_cell(Frame& frame);
    size_t min_numeric_type_ = 0;
    size_t max_numeric_type_ = 0;
    std::unordered_map<std::string, size_t> type_list_;
};

template <> const FormatMetadata& format_metadata<LAMMPSTrajectoryFormat>();

} // namespace chemfiles

#endif
