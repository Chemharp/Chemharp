// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_FORMAT_MOLFILE_HPP
#define CHEMFILES_FORMAT_MOLFILE_HPP

#include <cstddef>

#include <string>
#include <vector>

#include "chemfiles/File.hpp"
#include "chemfiles/Frame.hpp"
#include "chemfiles/Format.hpp"
#include "chemfiles/Topology.hpp"   // IWYU pragma: keep
#include "chemfiles/external/optional.hpp"

extern "C" {
    #include "vmdplugin.h"
    #include "molfile_plugin.h"
}

namespace chemfiles {
class FormatMetadata;

/// List all the VMD molfile plugins enabled. For more documentation about VMD
/// molfile plugins, please see:
/// http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/
enum MolfileFormat {
    TRJ,                ///< Gromacs .trj file format
    PSF,                ///< PSF topology files
    MOLDEN,             ///< Molden file format
};

/// A thin wrapper around the vmd plugin functions
template <MolfileFormat F>
struct MolfilePluginData final {
    /// Function to initialize the plugin
    int init();
    /// Function to register the plugin.
    int registration(void* data, vmdplugin_register_cb callback);
    /// Function to unload the plugin
    int fini();

    /// Plugin format
    std::string format() const;
    /// Plugin name
    std::string plugin_name() const;
    /// Plugin reader, when a given plugin manages different formats
    std::string reader() const;
    /// Does this plugin have velocity data
    bool have_velocities() const;
};

/// Use of VMD Molfile plugins as format reader. This class is templated by a
/// value in the `MolfileFormat` enum.
template <MolfileFormat F>
class Molfile final: public Format {
public:
    Molfile(std::string path, File::Mode mode, File::Compression compression);
    ~Molfile() override;

    void read(Frame& frame) override;
    void read_step(size_t step, Frame& frame) override;
    size_t nsteps() override;
private:
    /// Convert a molfile timestep to a chemfiles frame
    void molfile_to_frame(const molfile_timestep_t& timestep, Frame& frame);
    /// Read topological information in the current file, if any.
    void read_topology();

    int read_next_timestep(molfile_timestep_t* timestep);

    /// Path of the underlying file
    std::string path_;
    /// VMD plugin data
    MolfilePluginData<F> plugin_data_;
    /// VMD molfile plugin
    molfile_plugin_t* plugin_handle_;
    /// Data associated with the plugin, to be used by all the plugins callback
    void* data_;
    /// The number of atoms in this file
    int natoms_;
    /// Store optional topological information
    optional<Topology> topology_;
    /// Store pre-read steps to implement read_step
    std::vector<Frame> frames_;
};

template<> const FormatMetadata& format_metadata<Molfile<TRJ>>();
template<> const FormatMetadata& format_metadata<Molfile<PSF>>();
template<> const FormatMetadata& format_metadata<Molfile<MOLDEN>>();

} // namespace chemfiles

#endif
