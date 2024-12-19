// From vmdconio.h
#include <utility>
#define VMDCON_WARN      2
#define VMDCON_ERROR     3

extern "C" {
    #include <vmdplugin.h>
    #include <molfile_plugin.h>
}

#include <cassert>
#include <cstddef>
#include <cstdint>

#include <array>
#include <functional>
#include <string>
#include <vector>
#include <unordered_map>

#include "chemfiles/types.hpp"
#include "chemfiles/warnings.hpp"
#include "chemfiles/error_fmt.hpp"
#include "chemfiles/external/span.hpp"
#include "chemfiles/external/optional.hpp"

#include "chemfiles/File.hpp"
#include "chemfiles/Atom.hpp"
#include "chemfiles/Frame.hpp"
#include "chemfiles/Residue.hpp"
#include "chemfiles/Topology.hpp"
#include "chemfiles/FormatMetadata.hpp"

#include "chemfiles/formats/Molfile.hpp"

using namespace chemfiles;

/******************************************************************************/
#define PLUGINS_DATA(FORMAT, PLUGIN, READER, VELOCITIES)                       \
    extern "C" int PLUGIN##_register(void*, vmdplugin_register_cb); /*NOLINT*/ \
    extern "C" int PLUGIN##_fini(void); /*NOLINT*/                             \
    extern "C" int PLUGIN##_init(void); /*NOLINT*/                             \
    template <> struct MolfilePluginData<FORMAT> {                             \
        int init() { return PLUGIN##_init(); }                                 \
        int registration(void* data, vmdplugin_register_cb callback) {         \
            return PLUGIN##_register(data, callback);                          \
        }                                                                      \
        int fini() { return PLUGIN##_fini(); }                                 \
        std::string format() const {return #FORMAT;}                           \
        std::string plugin_name() const {return #PLUGIN;}                      \
        std::string reader() const {return #READER;}                           \
        bool have_velocities() const {return VELOCITIES;}                     \
    }

namespace chemfiles {
    PLUGINS_DATA(TRJ,     gromacsplugin, trj,    false);
    PLUGINS_DATA(PSF,     psfplugin,     psf,    false);
    PLUGINS_DATA(MOLDEN,  moldenplugin,  molden, false);
}

#undef PLUGINS_FUNCTIONS
/******************************************************************************/

/// data identifying a residue in molfile plugins
struct residue_info_t {
    int id;
    std::string name;
    std::string segid;
    std::string chain;

    bool operator==(const residue_info_t& other) const {
        return (
            this->id == other.id &&
            this->name == other.name &&
            this->segid == other.segid &&
            this->chain == other.chain
        );
    }
};

inline void hash_combine(size_t& /*hash*/) { }

template <typename T, typename... Tail>
inline void hash_combine(size_t& seed, const T& v, Tail... tail) {
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    hash_combine(seed, tail...);
}

namespace std {
    template<> struct hash<residue_info_t> {
        size_t operator()(const residue_info_t &info) const {
            size_t result = 0;
            hash_combine(result, info.id, info.name, info.segid, info.chain);
            return result;
        }
    };
}

template <MolfileFormat F> static int register_plugin(void* user_data, vmdplugin_t* vmd_plugin) {
    auto* user_plugin = static_cast<molfile_plugin_t**>(user_data);
    assert(std::string(MOLFILE_PLUGIN_TYPE) == std::string(vmd_plugin->type));

    auto* plugin = reinterpret_cast<molfile_plugin_t*>(vmd_plugin);
    // When this callback is called multiple times with more than one plugin,
    // only register the one whe want
    if (MolfilePluginData<F>().reader() == plugin->name) {
        *user_plugin = plugin;
    }

    return VMDPLUGIN_SUCCESS;
}

static int molfiles_to_chemfiles_warning(int level, const char* message) {
    if (level == VMDCON_ERROR || level == VMDCON_WARN) {
        send_warning(message);
    }
    return 0;
}

/******************************************************************************/

template <MolfileFormat F>
Molfile<F>::Molfile(std::string path, File::Mode mode, File::Compression compression)
    : path_(std::move(path)), plugin_handle_(nullptr), data_(nullptr), natoms_(0) {
    if (mode != File::READ) {
        throw format_error(
            "molfiles based format {} is only available in read mode", plugin_data_.format()
        );
    }

    if (compression != File::DEFAULT) {
        throw format_error(
            "molfiles based format {} do not support compression", plugin_data_.format()
        );
    }

    if (plugin_data_.init()) {
        throw format_error(
            "could not initialize the {} plugin", plugin_data_.format()
        );
    }

    if (plugin_data_.registration(&plugin_handle_, register_plugin<F>)) {
        throw format_error(
            "could not register the {} plugin", plugin_data_.format()
        );
    }

    // Check the ABI version of the plugin
    assert(plugin_handle_->abiversion == vmdplugin_ABIVERSION);
    // Redirect console output to chemfiles warnings
    plugin_handle_->cons_fputs = molfiles_to_chemfiles_warning;

    // Check that needed functions are here
    if (plugin_handle_->open_file_read == nullptr ||
        (
            plugin_handle_->read_next_timestep == nullptr &&
            plugin_handle_->read_timestep == nullptr &&
            plugin_handle_->read_structure == nullptr
        ) || plugin_handle_->close_file_read == nullptr ) {
        throw format_error(
            "the {} plugin does not have read capacities", plugin_data_.format()
        );
    }

    data_ = plugin_handle_->open_file_read(
        path_.c_str(), plugin_handle_->name, &natoms_
    );

    if (data_ == nullptr) {
        throw format_error(
            "could not open the file at '{}' with {} plugin", path_, plugin_data_.format()
        );
    }

    read_topology();
}

template <MolfileFormat F> Molfile<F>::~Molfile() {
    if (data_ != nullptr) {
        plugin_handle_->close_file_read(data_);
    }
    plugin_data_.fini();
}

template <MolfileFormat F> int Molfile<F>::read_next_timestep(molfile_timestep_t* timestep) {
    if (plugin_handle_->read_next_timestep != nullptr) {
        // This function is provided by classical molecular simulation format.
        return plugin_handle_->read_next_timestep(data_, natoms_, timestep);
    } else if (plugin_handle_->read_timestep != nullptr) {
        // This function is provided by quantum molecular simulation format.
        return plugin_handle_->read_timestep(
            data_, natoms_, timestep, nullptr, nullptr
        );
    } else if (plugin_handle_->read_structure != nullptr) {
        // topology only format, nothing to do
        return MOLFILE_SUCCESS;
    } else {
        throw format_error(
            "read_next_timestep, read_timestep and read_structure are missing "
            "in this plugin. This is a bug"
        );
    }
}

template <MolfileFormat F> void Molfile<F>::read(Frame& frame) {
    std::vector<float> coords(3 * static_cast<size_t>(natoms_));
    std::vector<float> velocities(0);

    molfile_timestep_t timestep{nullptr, nullptr, 0, 0, 0, 90, 90, 90, 0};
    timestep.coords = coords.data();
    if (plugin_data_.have_velocities()) {
        velocities.resize(3 * static_cast<size_t>(natoms_));
        timestep.velocities = velocities.data();
    }

    int status = read_next_timestep(&timestep);
    if (status != MOLFILE_SUCCESS) {
        throw format_error(
            "error while reading the file at '{}' with {} plugin",
            path_, plugin_data_.format()
        );
    }

    if (topology_) {
        frame.resize(topology_->size());
        frame.set_topology(*topology_);
    }
    molfile_to_frame(timestep, frame);

    frames_.emplace_back(frame.clone());
}

template <MolfileFormat F> void Molfile<F>::read_step(size_t step, Frame& frame) {
    while (step >= frames_.size()) {
        Frame new_frame;
        this->read(new_frame);
    }
    frame = frames_.at(step).clone();
}

template <MolfileFormat F> size_t Molfile<F>::nsteps() {
    if (plugin_handle_->read_next_timestep == nullptr) {
        // FIXME: this is hacky, but the molden plugin does not respect a NULL
        // argument for molfile_timestep_t, so for now we are only able to read
        // a single step from all the QM format plugins.
        return 1;
    }
    size_t n = 0;
    int status = MOLFILE_SUCCESS;
    while (true) {
        status = read_next_timestep(nullptr);
        if (status == MOLFILE_SUCCESS) {
            n++;
        } else {
            break;
        }
    }
    // We need to close and re-open the file
    plugin_handle_->close_file_read(data_);
    int unused = 0;
    data_ = plugin_handle_->open_file_read(path_.c_str(), plugin_handle_->name, &unused);
    read_topology();

    return n;
}

template <MolfileFormat F>
void Molfile<F>::molfile_to_frame(const molfile_timestep_t& timestep, Frame& frame) {
    auto lengths = Vector3D(
        static_cast<double>(timestep.A),
        static_cast<double>(timestep.B),
        static_cast<double>(timestep.C)
    );
    auto angles = Vector3D(
        static_cast<double>(timestep.alpha),
        static_cast<double>(timestep.beta),
        static_cast<double>(timestep.gamma)
    );
    frame.set_cell({lengths, angles});

    frame.resize(static_cast<size_t>(natoms_));
    auto positions = frame.positions();
    for (size_t i = 0; i < static_cast<size_t>(natoms_); i++) {
        positions[i][0] = static_cast<double>(timestep.coords[3 * i + 0]);
        positions[i][1] = static_cast<double>(timestep.coords[3 * i + 1]);
        positions[i][2] = static_cast<double>(timestep.coords[3 * i + 2]);
    }

    if (plugin_data_.have_velocities()) {
        frame.add_velocities();
        auto velocities = frame.velocities();
        for (size_t i = 0; i < static_cast<size_t>(natoms_); i++) {
            (*velocities)[i][0] = static_cast<double>(timestep.velocities[3 * i + 0]);
            (*velocities)[i][1] = static_cast<double>(timestep.velocities[3 * i + 1]);
            (*velocities)[i][2] = static_cast<double>(timestep.velocities[3 * i + 2]);
        }
    }
}

template <MolfileFormat F> void Molfile<F>::read_topology() {
    if (plugin_handle_->read_structure == nullptr) {
        return;
    }

    std::vector<molfile_atom_t> atoms(static_cast<size_t>(natoms_));
    int optflags = 0;
    int status = plugin_handle_->read_structure(data_, &optflags, atoms.data());
    if (status != MOLFILE_SUCCESS) {
        throw format_error(
            "could not read the molecule structure with {} plugin",
            plugin_data_.format()
        );
    }

    topology_ = Topology();

    auto residues = std::unordered_map<residue_info_t, Residue>();
    size_t atom_id = 0;
    for (auto& molfile_atom : atoms) {
        Atom atom(molfile_atom.name, molfile_atom.type);
        if ((optflags & MOLFILE_MASS) != 0) {
            atom.set_mass(static_cast<double>(molfile_atom.mass));
        }
        if ((optflags & MOLFILE_CHARGE) != 0) {
            atom.set_charge(static_cast<double>(molfile_atom.charge));
        }

        topology_->add_atom(std::move(atom));

        if (molfile_atom.resname != std::string("")) {
            auto resid = static_cast<int64_t>(molfile_atom.resid);
            auto residue = Residue(molfile_atom.resname, resid);
            residue.set("segname", molfile_atom.segid);
            residue.set("chainname", molfile_atom.chain);
            residue.set("chainid", molfile_atom.chain);
            auto info = residue_info_t {
                molfile_atom.resid,
                molfile_atom.resname,
                molfile_atom.segid,
                molfile_atom.chain,
            };
            auto inserted = residues.insert({std::move(info), std::move(residue)});
            inserted.first->second.add_atom(atom_id);
        }
        atom_id++;
    }

    for (auto residue: residues) {
        topology_->add_residue(std::move(residue.second));
    }

    if (plugin_handle_->read_bonds == nullptr) {
        return;
    }

    int nbonds = 0;
    int* from = nullptr;
    int* to = nullptr;

    int dummy = 0;
    float* _float = nullptr;
    int* _int = nullptr;
    char** _char = nullptr;

    status = plugin_handle_->read_bonds(
        data_, &nbonds, &from, &to, &_float, &_int, &dummy, &_char
    );
    if (status != MOLFILE_SUCCESS) {
        throw format_error(
            "could not read bonds with {} plugin", plugin_data_.format()
        );
    }

    for (size_t i = 0; i < static_cast<size_t>(nbonds); i++) {
        // Indexes are 1-based in Molfile
        topology_->add_bond(static_cast<size_t>(from[i] - 1),
                            static_cast<size_t>(to[i]) - 1);
    }
}

// Instantiate all the templates
template class chemfiles::Molfile<TRJ>;
template class chemfiles::Molfile<PSF>;
template class chemfiles::Molfile<MOLDEN>;

template<> const FormatMetadata& chemfiles::format_metadata<Molfile<TRJ>>() {
    static FormatMetadata metadata;
    metadata.name = "TRJ";
    metadata.extension = ".trj";
    metadata.description = "GROMACS .trj binary format";
    metadata.reference = "http://manual.gromacs.org/archive/5.0.7/online/trj.html";

    metadata.read = true;
    metadata.write = false;
    metadata.memory = false;

    metadata.positions = true;
    metadata.velocities = false;
    metadata.unit_cell = false;
    metadata.atoms = false;
    metadata.bonds = false;
    metadata.residues = false;
    return metadata;
}

template<> const FormatMetadata& chemfiles::format_metadata<Molfile<PSF>>() {
    static FormatMetadata metadata;
    metadata.name = "PSF";
    metadata.extension = ".psf";
    metadata.description = "Protein Structure File text format";
    metadata.reference = "https://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node23.html";

    metadata.read = true;
    metadata.write = false;
    metadata.memory = false;

    metadata.positions = false;
    metadata.velocities = false;
    metadata.unit_cell = false;
    metadata.atoms = true;
    metadata.bonds = true;
    // FIXME: the molfile plugin does not read residue information, we should
    // add it when re-implementing a PSF reader.
    metadata.residues = false;
    return metadata;
}

template<> const FormatMetadata& chemfiles::format_metadata<Molfile<MOLDEN>>() {
    static FormatMetadata metadata;
    metadata.name = "Molden";
    metadata.extension = ".molden";
    metadata.description = "Molden text format";
    metadata.reference = "https://web.archive.org/web/20200531080022/http://cheminf.cmbi.ru.nl/molden/molden_format.html";

    metadata.read = true;
    metadata.write = false;
    metadata.memory = false;

    metadata.positions = true;
    metadata.velocities = false;
    metadata.unit_cell = false;
    metadata.atoms = true;
    metadata.bonds = false;
    metadata.residues = false;
    return metadata;
}
