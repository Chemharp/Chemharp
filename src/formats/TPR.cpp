// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cassert>
#include <cstddef>
#include <cstdint>

#include <algorithm>
#include <array>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "chemfiles/error_fmt.hpp"
#include "chemfiles/external/optional.hpp"
#include "chemfiles/external/span.hpp"
#include "chemfiles/periodic_table.hpp"
#include "chemfiles/warnings.hpp"

#include "chemfiles/Atom.hpp"
#include "chemfiles/File.hpp"
#include "chemfiles/FormatMetadata.hpp"
#include "chemfiles/Frame.hpp"
#include "chemfiles/Residue.hpp"

#include "chemfiles/files/XDRFile.hpp"
#include "chemfiles/formats/TPR.hpp"

// see <GMX>/src/gromacs/fileio/tpxio.cpp
#define TPR_TAG_RELEASE "release"

// GROMACS explains: Enum of values that describe the contents of a TPR file
// whose format matches a version number.
// see `tpxv` in <GMX>/src/gromacs/fileio/tpxio.cpp
class TPRVersion {
  public:
    enum TV : int {
        ComputationalElectrophysiology = 96,
        Use64BitRandomSeed,
        RestrictedBendingAndCombinedAngleTorsionPotentials,
        InteractiveMolecularDynamics,
        RemoveObsoleteParameters1,
        PullCoordTypeGeom,
        PullGeomDirRel,
        IntermolecularBondeds,
        CompElWithSwapLayerOffset,
        CompElPolyatomicIonsAndMultipleIonTypes,
        RemoveAdress,
        PullCoordNGroup,
        RemoveTwinRange,
        ReplacePullPrintCOM12,
        PullExternalPotential,
        GenericParamsForElectricField,
        AcceleratedWeightHistogram,
        RemoveImplicitSolvation,
        PullPrevStepCOMAsReference,
        MimicQMMM,
        PullAverage,
        GenericInternalParameters,
        VSite2FD,
        AddSizeField,
        StoreNonBondedInteractionExclusionGroup,
        VSite1,
        MTS,
        RemovedConstantAcceleration,
        TransformationPullCoord,
        SoftcoreGapsys,
        ReaddedConstantAcceleration,
        RemoveTholeRfac,
        RemoveAtomtypes,
        EnsembleTemperature,
        Count // This number is for the total number of versions
    };
};

// GROMACS explains:
// Version number of the file format written to run input
// files by this version of the code.
//
// The TPR version increases whenever the file format in the main
// development branch changes, due to an extension of the `TPRVersion` enum above.
// Backward compatibility for reading old run input files is maintained
// by checking this version number against that of the file and then using
// the correct code path.
static const int TPR_VERSION = TPRVersion::Count - 1;

// GROMACS explains:
// Value of current TPR generation to keep track of incompatible changes for older TPR versions.
//
// The generation should be incremented when editing the TOPOLOGY or HEADER of the TPR format.
// In particular, updating `FUNCTION_TYPE_UPDATES` or changing the fields of the TPR header
// often triggers such needs.
//
// This way GROMACS can maintain forward compatibility for all analysis tools and/or external
// programs that only need to know the atom/residue names, charges, and bond connectivity.
//
// For chemfiles, this means it is possible to read TPR files from future versions as long as they
// are the same generation as specified here.
// The generation is tracked by GROMACS with an enum.
// see `TpxGeneration` in in <GMX>/src/gromacs/fileio/tpxio.cpp
static const int TPR_GENERATION = 28;

// This number should be the most recent backwards incompatible version.
// I.e., if this number is 9, we cannot read TPR version 9 with this code.
static const int TPR_INCOMPATIBLE_VERSION = 57; // GMX4.0 has version 58

// see <GMX>/api/legacy/include/gromacs/topology/ifunc.h
class FunctionType {
  public:
    enum FT : size_t {
        BONDS,
        G96BONDS,
        MORSE,
        CUBICBONDS,
        CONNBONDS,
        HARMONIC,
        FENEBONDS,
        TABBONDS,
        TABBONDSNC,
        RESTRBONDS,
        ANGLES,
        G96ANGLES,
        RESTRANGLES,
        LINEAR_ANGLES,
        CROSS_BOND_BONDS,
        CROSS_BOND_ANGLES,
        UREY_BRADLEY,
        QUARTIC_ANGLES,
        TABANGLES,
        PDIHS,
        RBDIHS,
        RESTRDIHS,
        CBTDIHS,
        FOURDIHS,
        IDIHS,
        PIDIHS,
        TABDIHS,
        CMAP,
        GB12_NOLONGERUSED,
        GB13_NOLONGERUSED,
        GB14_NOLONGERUSED,
        GBPOL_NOLONGERUSED,
        NPSOLVATION_NOLONGERUSED,
        LJ14,
        COUL14,
        LJC14_Q,
        LJC_PAIRS_NB,
        LJ,
        BHAM,
        LJ_LR_NOLONGERUSED,
        BHAM_LR_NOLONGERUSED,
        DISPCORR,
        COUL_SR,
        COUL_LR_NOLONGERUSED,
        RF_EXCL,
        COUL_RECIP,
        LJ_RECIP,
        DPD,
        POLARIZATION,
        WATER_POL,
        THOLE_POL,
        ANHARM_POL,
        POSRES,
        FBPOSRES,
        DISRES,
        DISRESVIOL,
        ORIRES,
        ORIRESDEV,
        ANGRES,
        ANGRESZ,
        DIHRES,
        DIHRESVIOL,
        CONSTR,
        CONSTRNC,
        SETTLE,
        VSITE1,
        VSITE2,
        VSITE2FD,
        VSITE3,
        VSITE3FD,
        VSITE3FAD,
        VSITE3OUT,
        VSITE4FD,
        VSITE4FDN,
        VSITEN,
        COM_PULL,
        DENSITYFITTING,
        EQM,
        EPOT,
        EKIN,
        ETOT,
        ECONSERVED,
        TEMP,
        VTEMP_NOLONGERUSED,
        PDISPCORR,
        PRES,
        DVDL_CONSTR,
        DVDL,
        DKDL,
        DVDL_COUL,
        DVDL_VDW,
        DVDL_BONDED,
        DVDL_RESTRAINT,
        DVDL_TEMPERATURE,
        Count
    };
    FunctionType() = default;
    constexpr FunctionType(FT v) : value(v) {}
    constexpr FunctionType(size_t v) : value(static_cast<FT>(v)) {}
    FunctionType(int v) : value(static_cast<FT>(v)) { assert(v >= 0); }
    operator int() const { return static_cast<int>(value); }

  private:
    FT value;
};

// This number is for the total number of energies / interaction function types
static const int NRE = FunctionType::Count;

// Set of function types which are considered as bonds.
// Update when new function types are introduced.
const std::vector<FunctionType> BOND_TYPES = {
    FunctionType::BONDS,      FunctionType::G96BONDS,   FunctionType::MORSE,
    FunctionType::CUBICBONDS, FunctionType::CONNBONDS,  FunctionType::HARMONIC,
    FunctionType::FENEBONDS,  FunctionType::RESTRBONDS, FunctionType::CONSTR,
    FunctionType::CONSTRNC,   FunctionType::TABBONDS,   FunctionType::TABBONDSNC,
};

// All force filed parameters of the system.
// see `gmx_ffparams_t` in <GMX>/api/legacy/include/gromacs/topology/forcefieldparameters.h
struct FFParams {
    // Number of non-bonded atom types
    int natom_types;
    // Function type of an interaction
    std::vector<int32_t> function_types;
};

// Maximum allowed number of parameters.
// see <GMX>/api/legacy/include/gromacs/topology/ifunc.h
const int NR_RBDIHS = 6;
const int NR_CBTDIHS = 6;

// see `t_ftupd` in <GMX>/src/gromacs/fileio/tpxio.cpp
struct FunctionTypeUpdate {
    // File version number in which the function type first appeared
    int file_version;
    FunctionType function_type;
};

// The following Table is used maintain TRR compatibility  when function types
// are added. Necessary to support reading of old TPR file versions.
// see `ftupd` in <GMX>/src/gromacs/fileio/tpxio.cpp
static const FunctionTypeUpdate FUNCTION_TYPE_UPDATES[24]{
    {70, FunctionType::RESTRBONDS},
    {TPRVersion::RestrictedBendingAndCombinedAngleTorsionPotentials, FunctionType::RESTRANGLES},
    {76, FunctionType::LINEAR_ANGLES},
    {TPRVersion::RestrictedBendingAndCombinedAngleTorsionPotentials, FunctionType::RESTRDIHS},
    {TPRVersion::RestrictedBendingAndCombinedAngleTorsionPotentials, FunctionType::CBTDIHS},
    {65, FunctionType::CMAP},
    {60, FunctionType::GB12_NOLONGERUSED},
    {61, FunctionType::GB13_NOLONGERUSED},
    {61, FunctionType::GB14_NOLONGERUSED},
    {72, FunctionType::GBPOL_NOLONGERUSED},
    {72, FunctionType::NPSOLVATION_NOLONGERUSED},
    {93, FunctionType::LJ_RECIP},
    {76, FunctionType::ANHARM_POL},
    {90, FunctionType::FBPOSRES},
    {TPRVersion::VSite1, FunctionType::VSITE1},
    {TPRVersion::VSite2FD, FunctionType::VSITE2FD},
    {TPRVersion::GenericInternalParameters, FunctionType::DENSITYFITTING},
    {69, FunctionType::VTEMP_NOLONGERUSED},
    {66, FunctionType::PDISPCORR},
    {79, FunctionType::DVDL_COUL},
    {79, FunctionType::DVDL_VDW},
    {79, FunctionType::DVDL_BONDED},
    {79, FunctionType::DVDL_RESTRAINT},
    {79, FunctionType::DVDL_TEMPERATURE},
};

// see `t_interaction_function` in <GMX>/api/legacy/include/gromacs/topology/ifunc.h
struct FunctionTypeInfo {
    std::string name;
    size_t natoms;
};

// This table contains a human-readable name for a function type
// and the number of atoms needed per function.
// see `interaction_function` in <GMX>/src/gromacs/topology/ifunc.cpp
const FunctionTypeInfo FUNCTION_TYPE_INFOS[NRE]{
    {"Bond", 2},
    {"G96Bond", 2},
    {"Morse", 2},
    {"Cubic Bonds", 2},
    {"Connect Bonds", 2},
    {"Harmonic Pot.", 2},
    {"FENE Bonds", 2},
    {"Tab. Bonds", 2},
    {"Tab. Bonds NC", 2},
    {"Restraint Pot.", 2},
    {"Angle", 3},
    {"G96Angle", 3},
    {"Restricted Angles", 3},
    {"Lin. Angle", 3},
    {"Bond-Cross", 3},
    {"BA-Cross", 3},
    {"U-B", 3},
    {"Quartic Angles", 3},
    {"Tab. Angles", 3},
    {"Proper Dih.", 4},
    {"Ryckaert-Bell.", 4},
    {"Restricted Dih.", 4},
    {"CBT Dih.", 4},
    {"Fourier Dih.", 4},
    {"Improper Dih.", 4},
    {"Periodic Improper Dih.", 4},
    {"Tab. Dih.", 4},
    {"CMAP Dih.", 5},
    {"GB 1-2 Pol. ", 0},
    {"GB 1-3 Pol. ", 0},
    {"GB 1-4 Pol. ", 0},
    {"GB Polarization ", 0},
    {"Nonpolar Sol. ", 0},
    {"LJ-14", 2},
    {"Coulomb-14", 0},
    {"LJC-14 q", 2},
    {"LJC Pairs NB", 2},
    {"LJ (SR)", 2},
    {"Buck.ham (SR)", 2},
    {"LJ ", 0},
    {"B.ham ", 0},
    {"Disper. corr.", 0},
    {"Coulomb (SR)", 0},
    {"Coul ", 0},
    {"RF excl.", 0},
    {"Coul. recip.", 0},
    {"LJ recip.", 0},
    {"DPD", 0},
    {"Polarization", 2},
    {"Water Pol.", 5},
    {"Thole Pol.", 4},
    {"Anharm. Pol.", 2},
    {"Position Rest.", 1},
    {"Flat-bottom posres", 1},
    {"Dis. Rest.", 2},
    {"D.R.Viol. (nm)", 0},
    {"Orient. Rest.", 2},
    {"Ori. R. RMSD", 0},
    {"Angle Rest.", 4},
    {"Angle Rest. Z", 2},
    {"Dih. Rest.", 4},
    {"Dih. Rest. Viol.", 0},
    {"Constraint", 2},
    {"Constr. No Conn.", 2},
    {"Settle", 3},
    {"Virtual site 1", 2},
    {"Virtual site 2", 3},
    {"Virtual site 2fd", 3},
    {"Virtual site 3", 4},
    {"Virtual site 3fd", 4},
    {"Virtual site 3fad", 4},
    {"Virtual site 3out", 4},
    {"Virtual site 4fd", 5},
    {"Virtual site 4fdn", 5},
    {"Virtual site N", 2},
    {"COM Pull En.", 0},
    {"Density fitting", 0},
    {"Quantum En.", 0},
    {"Potential", 0},
    {"Kinetic En.", 0},
    {"Total Energy", 0},
    {"Conserved En.", 0},
    {"Temperature", 0},
    {"Vir. Temp. ", 0},
    {"Pres. DC", 0},
    {"Pressure", 0},
    {"dH/dl constr.", 0},
    {"dVremain/dl", 0},
    {"dEkin/dl", 0},
    {"dVcoul/dl", 0},
    {"dVvdw/dl", 0},
    {"dVbonded/dl", 0},
    {"dVrestraint/dl", 0},
    {"dVtemperature/dl", 0},
};

// see `t_resinfo` in <GMX>/api/legacy/include/gromacs/topology/atoms.h
struct ResidueInfo {
    std::string name;
    int64_t residue_number;
    uint8_t insertion_code;
};

// see `t_atom` in <GMX>/api/legacy/include/gromacs/topology/atoms.h
struct AtomProperties {
    double mass;
    double charge;
    // Index into `residue_infos` in `Atoms`
    size_t residue_idx;
    chemfiles::optional<std::string> element_name = chemfiles::nullopt;
};

// see `t_atoms` in <GMX>/api/legacy/include/gromacs/topology/atoms.h
struct Atoms {
    std::vector<AtomProperties> atom_properties;
    std::vector<std::string> atom_names;
    std::vector<std::string> atom_types;
    std::vector<ResidueInfo> residue_infos;

    size_t size() const {
        size_t natoms = atom_properties.size();
        assert(natoms == atom_names.size());
        assert(natoms == atom_types.size());
        return natoms;
    }
};

// All the interactions for a given function type. This is stored as multiple
// interaction tuples which refer to a set of function parameters (skipped)
// and the atom indices.
// see <GMX>/api/legacy/include/gromacs/topology/idef.h
struct InteractionList {
    FunctionType function_type;
    // An interaction tuple contains the type of the interaction
    // (i.e. index to function parameters) first
    // followed by the number of interacting atoms.
    // This number depends on the function type (e.g. bond defines two atoms).
    std::vector<size_t> interaction_tuples;
    bool settle_done = false;

    bool empty() const { return interaction_tuples.empty(); }

    size_t size() const {
        // Number of entries per interaction (type_index + N * atom_index)
        size_t nentries = FUNCTION_TYPE_INFOS[function_type].natoms + 1;
        size_t size = interaction_tuples.size() / nentries;
        assert(size * nentries == interaction_tuples.size());
        return size;
    }

    chemfiles::span<const size_t> operator[](size_t i) const {
        size_t natoms = FUNCTION_TYPE_INFOS[function_type].natoms;
        // Stored as (type, atom_1, atom_2, ..., atom_N)
        // Skip first `i` tuples and the next type
        size_t begin = (natoms + 1) * i + 1;
        return chemfiles::span<const size_t>(interaction_tuples.data() + begin, natoms);
    }

    void add_settle_atoms() {
        assert(function_type == FunctionType::SETTLE);
        assert(!settle_done); // This can be called only once
        settle_done = true;
        // GROMACS: Settle used to only store the first atom: add the other two
        size_t size = interaction_tuples.size();
        interaction_tuples.resize(2 * size);
        auto& iatoms = interaction_tuples;
        for (int i = static_cast<int>(size) / 2 - 1; i >= 0; --i) {
            size_t j = static_cast<size_t>(i);
            iatoms[4 * j + 0] = iatoms[2 * j + 0];
            iatoms[4 * j + 1] = iatoms[2 * j + 1];
            iatoms[4 * j + 2] = iatoms[2 * j + 1] + 1;
            iatoms[4 * j + 3] = iatoms[2 * j + 1] + 2;
        }
    }
};

// GROMACS: List of interaction lists, one list for each interaction type.
// see <GMX>/api/legacy/include/gromacs/topology/idef.h
typedef std::array<chemfiles::optional<InteractionList>, NRE> InteractionLists;

// see `gmx_moltype_t` in <GMX>/api/legacy/include/gromacs/topology/topology.h
struct MoleculeType {
    std::string name;
    Atoms atoms;
    InteractionLists interaction_lists;
};

// `Count` in `SimulationAtomGroupType`
// see <GMX>/api/legacy/include/gromacs/topology/topology_enums.h
const size_t NR_GROUP_TYPES = 10;

using namespace chemfiles;

template <> const FormatMetadata& chemfiles::format_metadata<TPRFormat>() {
    static FormatMetadata metadata;
    metadata.name = "TPR";
    metadata.extension = ".tpr";
    metadata.description = "GROMACS TPR binary format";
    metadata.reference = "http://manual.gromacs.org/current/reference-manual/file-formats.html#tpr";

    metadata.read = true;
    metadata.write = false;
    metadata.memory = false;

    metadata.positions = true;
    metadata.velocities = true;
    metadata.unit_cell = true;
    metadata.atoms = true;
    metadata.bonds = true;
    metadata.residues = true;
    return metadata;
}

TPRFormat::TPRFormat(std::string path, File::Mode mode, File::Compression compression):
    file_(std::move(path), mode)
{
    if (compression != File::DEFAULT) {
        throw format_error("TPR format does not support compression");
    }
    if (mode != File::READ) {
        throw format_error("TPR format does not support write & append");
    }
    read_header();
}

size_t TPRFormat::nsteps() { return 1; }

void TPRFormat::read_step(size_t step, Frame& frame) {
    step_ = step;
    this->read(frame);
}

// Determine the size of all the interaction parameters for a function type.
// see `do_iparams` <GMX>/src/gromacs/fileio/tpxio.cpp
static size_t interaction_params_size(FunctionType ftype, size_t sizeof_real, int file_version) {
    switch (ftype) {
    case FunctionType::ANGLES:
    case FunctionType::G96ANGLES:
    case FunctionType::BONDS:
    case FunctionType::G96BONDS:
    case FunctionType::HARMONIC:
    case FunctionType::IDIHS:
        return 4 * sizeof_real;
    case FunctionType::RESTRANGLES:
        return 2 * sizeof_real;
    case FunctionType::LINEAR_ANGLES:
        return 4 * sizeof_real;
    case FunctionType::FENEBONDS:
        return 2 * sizeof_real;
    case FunctionType::RESTRBONDS:
        return 8 * sizeof_real;
    case FunctionType::TABBONDS:
    case FunctionType::TABBONDSNC:
    case FunctionType::TABANGLES:
    case FunctionType::TABDIHS:
        return 2 * sizeof_real + sizeof(int32_t);
    case FunctionType::CROSS_BOND_BONDS:
        return 3 * sizeof_real;
    case FunctionType::CROSS_BOND_ANGLES:
        return 4 * sizeof_real;
    case FunctionType::UREY_BRADLEY: {
        size_t size = 4 * sizeof_real;
        if (file_version >= 79) {
            size += 4 * sizeof_real;
        }
        return size;
    }
    case FunctionType::QUARTIC_ANGLES:
        return 6 * sizeof_real;
    case FunctionType::BHAM:
        return 3 * sizeof_real;
    case FunctionType::MORSE: {
        size_t size = 3 * sizeof_real;
        if (file_version >= 79) {
            size += 3 * sizeof_real;
        }
        return size;
    }
    case FunctionType::CUBICBONDS:
        return 3 * sizeof_real;
    case FunctionType::CONNBONDS:
        return 0;
    case FunctionType::POLARIZATION:
        return sizeof_real;
    case FunctionType::ANHARM_POL:
        return 3 * sizeof_real;
    case FunctionType::WATER_POL:
        return 6 * sizeof_real;
    case FunctionType::THOLE_POL: {
        size_t size = 3 * sizeof_real;
        if (file_version < TPRVersion::RemoveTholeRfac) {
            size += sizeof_real;
        }
        return size;
    }
    case FunctionType::LJ:
        return 2 * sizeof_real;
    case FunctionType::LJ14:
        return 4 * sizeof_real;
    case FunctionType::LJC14_Q:
        return 5 * sizeof_real;
    case FunctionType::LJC_PAIRS_NB:
        return 4 * sizeof_real;
    case FunctionType::PDIHS:
    case FunctionType::PIDIHS:
    case FunctionType::ANGRES:
    case FunctionType::ANGRESZ:
        return 4 * sizeof_real + sizeof(int32_t);
    case FunctionType::RESTRDIHS:
        return 2 * sizeof_real;
    case FunctionType::DISRES:
        return 2 * sizeof(int32_t) + 4 * sizeof_real;
    case FunctionType::ORIRES:
        return 3 * sizeof(int32_t) + 3 * sizeof_real;
    case FunctionType::DIHRES: {
        size_t size = 3 * sizeof_real;
        if (file_version < 82) {
            size += 2 * sizeof(int32_t);
        }
        if (file_version >= 82) {
            size += 3 * sizeof_real;
        }
        return size;
    }
    case FunctionType::POSRES:
        return 4 * 3 * sizeof_real;
    case FunctionType::FBPOSRES:
        return sizeof(int32_t) + 5 * sizeof_real;
    case FunctionType::CBTDIHS:
        return NR_CBTDIHS * sizeof_real;
    case FunctionType::RBDIHS:
    case FunctionType::FOURDIHS:
        return 2 * NR_RBDIHS * sizeof_real;
    case FunctionType::CONSTR:
    case FunctionType::CONSTRNC:
    case FunctionType::SETTLE:
        return 2 * sizeof_real;
    case FunctionType::VSITE1:
        return 0;
    case FunctionType::VSITE2:
    case FunctionType::VSITE2FD:
        return 1 * sizeof_real;
    case FunctionType::VSITE3:
    case FunctionType::VSITE3FD:
    case FunctionType::VSITE3FAD:
        return 2 * sizeof_real;
    case FunctionType::VSITE3OUT:
    case FunctionType::VSITE4FD:
    case FunctionType::VSITE4FDN:
        return 3 * sizeof_real;
    case FunctionType::VSITEN:
        return sizeof(int32_t) + sizeof_real;
    case FunctionType::GB12_NOLONGERUSED:
    case FunctionType::GB13_NOLONGERUSED:
    case FunctionType::GB14_NOLONGERUSED: {
        size_t size = 0;
        if (file_version < 68) {
            size += 4 * sizeof_real;
        }
        if (file_version < TPRVersion::RemoveImplicitSolvation) {
            size += 5 * sizeof_real;
        }
        return size;
    }
    case FunctionType::CMAP:
        return 2 * sizeof(int32_t);
    default:
        throw format_error(
            "unknown function type {} ({})",
            static_cast<int>(ftype),
            FUNCTION_TYPE_INFOS[ftype].name
        );
    }
}

// see `do_ilists` in <GMX>/src/gromacs/fileio/tpxio.cpp
static InteractionLists read_interaction_lists(XDRFile& file, int file_version) {
    InteractionLists interaction_lists;
    for (size_t i = 0; i < NRE; ++i) {
        const FunctionType function_type = static_cast<FunctionType>(i);
        bool is_old_interaction = false;
        for (const auto function_type_update : FUNCTION_TYPE_UPDATES) {
            // Compare the read file version to the update table
            if ((file_version < function_type_update.file_version) &&
                (function_type == function_type_update.function_type)) {
                is_old_interaction = true;
                break;
            }
        }
        if (!is_old_interaction) {
            InteractionList ilist;

            ilist.function_type = function_type;
            size_t ninteraction_data = file.read_single_size_as_i32();
            ilist.interaction_tuples.reserve(ninteraction_data);
            for (size_t j = 0; j < ninteraction_data; ++j) {
                size_t idx = file.read_single_size_as_i32();
                ilist.interaction_tuples.emplace_back(idx);
            }

            if (file_version < 78 && ilist.function_type == FunctionType::SETTLE &&
                !ilist.empty()) {
                ilist.add_settle_atoms();
            }

            interaction_lists[i] = ilist;
        }
    }
    return interaction_lists;
}

// Add connectivity elements i.e. bonds to the frame.
// Use the atom index offset to correct for molecule-internal numbering.
static void add_conectivity(Frame& frame, const InteractionLists& interaction_lists,
                            size_t atom_idx_offset = 0) {
    auto contains = [](const std::vector<FunctionType>& types_set,
                       FunctionType function_type) -> bool {
        auto cnt = std::count(types_set.begin(), types_set.end(), function_type);
        return cnt != 0;
    };
    for (const auto& ilist : interaction_lists) {
        if (!ilist) {
            continue;
        } else if (contains(BOND_TYPES, ilist.value().function_type)) {
            for (size_t i = 0; i < ilist.value().size(); ++i) {
                auto iatoms = ilist.value()[i];
                assert(iatoms.size() == 2);
                frame.add_bond(atom_idx_offset + iatoms[0], atom_idx_offset + iatoms[1]);
            }
        } else if (ilist.value().function_type == FunctionType::SETTLE) {
            for (size_t i = 0; i < ilist.value().size(); ++i) {
                auto iatoms = ilist.value()[i];
                assert(iatoms.size() == 3);
                frame.add_bond(atom_idx_offset + iatoms[0], atom_idx_offset + iatoms[1]);
                frame.add_bond(atom_idx_offset + iatoms[0], atom_idx_offset + iatoms[2]);
            }
        }
    }
}

// Read all the different interations in the system. As only the types of
// interactions are used, the interaction parameters are skipped.
// see `do_ffparams` in <GMX>/src/gromacs/fileio/tpxio.cpp
static FFParams read_force_field_parameters(XDRFile& file, size_t sizeof_real, int file_version) {
    FFParams ffparams;
    ffparams.natom_types = file.read_single_i32();
    const size_t num_function_types = file.read_single_size_as_i32();
    ffparams.function_types.resize(num_function_types);

    // Read all function types
    file.read_i32(ffparams.function_types);

    if (file_version >= 66) {
        // Skip `reppow`: Repulsion power $p$ for VdW: $C12*r^-p$
        file.skip(sizeof(double));
    }

    // Skip `fudgeQQ`: Scaling factor $f$ for Coulomb 1-4: $f*q1*q2$
    file.skip(sizeof_real);

    // GROMACS explains:
    // Check whether all these function types are supported by the code.
    // In practice the code is backwards compatible, which means that the
    // numbering may have to be altered from old numbering to new numbering.
    for (auto& function_type : ffparams.function_types) {
        for (const auto& function_type_update : FUNCTION_TYPE_UPDATES) {
            // Compare the read file version to the update table
            if ((file_version < function_type_update.file_version) &&
                (function_type >= function_type_update.function_type)) {
                ++function_type;
            }
        }
        // Skip the force field parameters of the function
        size_t params_size = interaction_params_size(static_cast<FunctionType>(function_type),
                                                     sizeof_real, file_version);
        file.skip(params_size);
    }
    return ffparams;
}

void TPRFormat::read(Frame& frame) {
    frame.resize(header_.natoms);

    // Now read the body of the TPR file
    // see `do_tpx_body` in <GMX>/src/gromacs/fileio/tpxio.cpp

    if (header_.has_box) {
        read_box(frame);
    }

    if (header_.ngroups_temperature_coupling > 0) {
        const size_t ngtc_size = header_.ngroups_temperature_coupling * header_.sizeof_real();
        if (header_.file_version < 69) {
            // Skip some legacy entries
            file_.skip(ngtc_size);
        }
        // GROMACS: These used to be the Berendsen tcoupl_lambda's
        file_.skip(ngtc_size);
    }

    if (header_.has_topology) {
        read_topology(frame);
    }

    read_coordinates(frame);

    // Now comes the simulation parameters in the input record
    // The parsing is done by `do_tpx_finalize()`
    // We stop parsing here

    step_++;
}

void TPRFormat::read_header() {
    const std::string version = file_.read_gmx_string();
    if (version.compare(0, 7, "VERSION") != 0) {
        throw format_error("unsupported file from a GROMACS version which is older than 2.0");
    }

    size_t precision = file_.read_single_size_as_i32();
    header_.use_double = (precision == sizeof(double));
    if (precision != sizeof(float) && precision != sizeof(double)) {
        throw format_error("invalid precision {}, expected {} or {}", precision, sizeof(float),
                           sizeof(double));
    }

    header_.file_version = file_.read_single_i32();

    // GROMACS explains:
    // This is for backward compatibility with development versions 77-79
    // where the tag was, mistakenly, placed before the generation,
    // which would cause a segv instead of a proper error message
    // when reading the topology only from tpx with <77 code.
    std::string fileTag;
    if (header_.file_version >= 77 && header_.file_version <= 79) {
        fileTag = file_.read_gmx_string();
    }

    header_.file_generation = file_.read_single_i32();

    if (header_.file_version >= 81) {
        fileTag = file_.read_gmx_string();
    }
    if (header_.file_version < 77 || header_.file_version == 80) {
        // GROMACS: Versions before 77 don't have the tag, set it to release.
        // Version 80 is not handled by the current GROMACS implementation
        // but MDAnalysis sets the tag to release as well for version 80.
        // Version 80 was used by both 5.0-dev and 4.6-dev.
        fileTag = TPR_TAG_RELEASE;
    }

    // GROMACS explains:
    // We only support reading tpx files with the same tag as the code
    // or tpx files with the release tag and with lower version number.
    if (fileTag != TPR_TAG_RELEASE && header_.file_version < TPR_VERSION) {
        throw format_error("TPR tag/version mismatch: reading file with version {}, tag '{}' with "
                           "program for version {}, tag '{}'",
                           header_.file_version, fileTag, TPR_VERSION, TPR_TAG_RELEASE);
    }

    if (header_.file_version > TPR_VERSION) {
        warning("TPR",
                "file version is from the future: implementation uses version {}, but got {}",
                TPR_VERSION, header_.file_version);
    }

    // Assume `TopOnlyOK` is true, i.e. we need only the topology and not the input record.
    // This allows reading of future versions of the same generation.
    if ((header_.file_version <= TPR_INCOMPATIBLE_VERSION) ||
        (header_.file_generation > TPR_GENERATION)) {
        throw format_error("unable to read version {} with version {} program",
                           header_.file_version, TPR_VERSION);
    }

    header_.natoms = file_.read_single_size_as_i32();
    header_.ngroups_temperature_coupling = file_.read_single_size_as_i32();

    if (header_.file_version < 62) {
        // Skip some legacy entries
        file_.skip(sizeof(int) + header_.sizeof_real());
    }
    if (header_.file_version >= 79) {
        // Skip current value of the alchemical state
        file_.read_single_i32();
    }
    if (!header_.use_double) {
        header_.lambda = static_cast<double>(file_.read_single_f32());
    } else {
        header_.lambda = file_.read_single_f64();
    }
    header_.has_input_record = read_gmx_bool();
    header_.has_topology = read_gmx_bool();
    header_.has_positions = read_gmx_bool();
    header_.has_velocities = read_gmx_bool();
    header_.has_forces = read_gmx_bool();
    header_.has_box = read_gmx_bool();

    if (header_.file_version >= TPRVersion::AddSizeField && header_.file_generation >= 27) {
        // Skip size of the TPR body in bytes
        file_.read_single_i64();
    }

    if (header_.file_generation > TPR_GENERATION && header_.has_input_record) {
        // Trying to read a file from the future with an input record.
        // At this point in time, it's unknown what will be in the record.
        header_.has_input_record = false;
    }

    if (header_.file_version >= TPRVersion::AddSizeField && header_.file_generation >= 27) {
        header_.body_convention = InMemory;
    } else {
        header_.body_convention = FileIOXdr;
    }
}

void TPRFormat::read_box(Frame& frame) {
    const auto box = file_.read_gmx_box(header_.use_double);
    frame.set_cell(box);

    if (header_.file_version >= 51) {
        // Relative box vectors characteristic of the box shape
        // Skip unused 3*3 real matrix
        file_.skip(header_.sizeof_real() * 9);
    }

    // Box velocities for Parrinello-Rahman barostat
    // Skip unused 3*3 real matrix
    file_.skip(header_.sizeof_real() * 9);

    if (header_.file_version < 56) {
        // Skip some legacy entries
        file_.skip(header_.sizeof_real() * 9);
    }
}

void TPRFormat::read_topology(Frame& frame) {
    // Most of the used functions are defined in <GMX>/src/gromacs/fileio/tpxio.cpp
    // which is not repeated everytime in the following section.

    const std::vector<std::string> symbol_table = read_symbol_table();

    frame.set("name", read_symbol_table_entry(symbol_table));

    const FFParams ffparams =
        read_force_field_parameters(file_, header_.sizeof_real(), header_.file_version);

    // Read the definitions of the different molecule types and
    // their atoms and residues.
    // see `do_moltype`
    auto read_molecule_types = [this, symbol_table]() -> std::vector<MoleculeType> {
        const size_t nmoltypes = file_.read_single_size_as_i32();
        std::vector<MoleculeType> molecule_types;
        molecule_types.reserve(nmoltypes);

        // Read the properties of single atom
        // see `do_atom`
        auto read_atom = [this]() -> AtomProperties {
            AtomProperties atom_prop;
            if (!header_.use_double) {
                // Float
                atom_prop.mass = static_cast<double>(file_.read_single_f32());
                atom_prop.charge = static_cast<double>(file_.read_single_f32());
            } else {
                // Double
                atom_prop.mass = file_.read_single_f64();
                atom_prop.charge = file_.read_single_f64();
            }
            // Skip mass and charge for Free Energy calculations
            file_.skip(2 * header_.sizeof_real());
            // Skip internal atom type
            if (header_.body_convention == FileIOXdr) {
                file_.skip(2 * sizeof(uint32_t));
            } else {
                assert(header_.body_convention == InMemory);
                file_.skip(2 * sizeof(uint16_t));
            }
            // Skip internal particle type
            file_.skip(sizeof(int32_t));
            atom_prop.residue_idx = file_.read_single_size_as_i32();
            const uint64_t atomnumber = static_cast<uint64_t>(file_.read_single_i32());
            for (const auto& kv : PERIODIC_TABLE) {
                // Determine the shortened element name from the atomnumber
                const auto& name = kv.first;
                const auto& atomic_data = kv.second;
                if (atomic_data.number && *atomic_data.number == atomnumber) {
                    atom_prop.element_name = name;
                }
            }
            return atom_prop;
        };

        // Read the definition of all the atoms in a molecule type
        // -> do_atoms
        auto read_atoms = [this, symbol_table, read_atom]() -> Atoms {
            Atoms atoms;
            const size_t natoms = file_.read_single_size_as_i32();
            const size_t nresidue_info = file_.read_single_size_as_i32();
            atoms.atom_properties.reserve(natoms);
            atoms.atom_names.reserve(natoms);
            atoms.atom_types.reserve(natoms);
            atoms.residue_infos.reserve(nresidue_info);
            for (size_t j = 0; j < natoms; ++j) {
                AtomProperties atom_prop = read_atom();
                atoms.atom_properties.emplace_back(atom_prop);
            }
            for (size_t j = 0; j < natoms; ++j) {
                auto name = read_symbol_table_entry(symbol_table);
                atoms.atom_names.emplace_back(name);
            }
            for (size_t j = 0; j < natoms; ++j) {
                auto type = read_symbol_table_entry(symbol_table);
                atoms.atom_types.emplace_back(type);
            }
            for (size_t j = 0; j < natoms; ++j) {
                // Skip the types for Free Energy calculations
                read_symbol_table_entry(symbol_table);
            }

            // Read the definition of residues
            // see `do_resinfo`
            for (size_t j = 0; j < nresidue_info; ++j) {
                std::string name = read_symbol_table_entry(symbol_table);
                uint8_t insertion_code = ' ';
                int64_t residue_number = static_cast<int64_t>(j) + 1;
                if (header_.file_version >= 63) {
                    residue_number = static_cast<int64_t>(file_.read_single_i32());
                    insertion_code = read_gmx_uchar();
                }
                ResidueInfo res_info = {name, residue_number, insertion_code};
                atoms.residue_infos.emplace_back(res_info);
            }

            return atoms;
        };

        for (size_t i = 0; i < nmoltypes; ++i) {
            MoleculeType moltype;
            moltype.name = read_symbol_table_entry(symbol_table);
            moltype.atoms = read_atoms();

            // Read the interaction lists
            moltype.interaction_lists = read_interaction_lists(file_, header_.file_version);

            molecule_types.emplace_back(moltype);

            // GROMACS: Skip the obsolete charge group index
            // see `do_block`
            const size_t nblocks = file_.read_single_size_as_i32();
            file_.skip((nblocks + 1) * sizeof(int32_t));

            // Skip the exclusions
            size_t nlists = file_.read_single_size_as_i32();
            size_t nelements = file_.read_single_size_as_i32();
            file_.skip((nlists + 1 + nelements) * sizeof(int32_t));
        }
        return molecule_types;
    };
    const std::vector<MoleculeType> molecule_types = read_molecule_types();

    // Now create the actual topology of the system. The molecules are listed
    // sequentially after each other. Multiple molecules of the same type in
    // one row are aggregated in molecule blocks.
    // see `do_molblock` but most of the code is chemfiles specific
    size_t global_atom_idx = 0; // Number of atoms in the previous molecules
    const size_t nmolblocks = file_.read_single_size_as_i32();
    for (size_t i = 0; i < nmolblocks; ++i) {
        // Index of the molecule type read previously
        size_t moltype_idx = file_.read_single_size_as_i32();
        const auto& moltype = molecule_types.at(moltype_idx);
        size_t nmolecules = file_.read_single_size_as_i32();
        const auto& atoms = moltype.atoms;
        for (size_t j = 0; j < nmolecules; ++j) {
            std::vector<chemfiles::Residue> residues_of_mol;
            residues_of_mol.reserve(atoms.residue_infos.size());
            for (const auto& res_info : atoms.residue_infos) {
                Residue residue = {res_info.name, res_info.residue_number};
                residue.set("insertion_code",
                            std::string(1, static_cast<char>(res_info.insertion_code)));
                residues_of_mol.emplace_back(residue);
            }
            for (size_t atom_idx = 0; atom_idx < atoms.size(); ++atom_idx) {
                auto& atom = frame[global_atom_idx + atom_idx];
                atom.set_name(atoms.atom_names[atom_idx]);
                atom.set("ff_type", atoms.atom_types[atom_idx]);
                const auto& props = atoms.atom_properties[atom_idx];
                if (props.element_name) {
                    atom.set_type(*props.element_name);
                }
                atom.set_mass(props.mass);
                atom.set_charge(props.charge);

                if (props.residue_idx < moltype.atoms.residue_infos.size()) {
                    residues_of_mol[props.residue_idx].add_atom(global_atom_idx + atom_idx);
                } else {
                    throw format_error(
                        "residue index out of bounds, there are {} residues, got index {}",
                        moltype.atoms.residue_infos.size(), props.residue_idx);
                }

                add_conectivity(frame, moltype.interaction_lists, global_atom_idx);
            }
            global_atom_idx += atoms.size();
            for (const auto& residue : residues_of_mol) {
                frame.add_residue(residue);
            }
        }

        // Skip number of atoms per molecule as this is already stored in the type
        file_.skip(sizeof(int32_t));
        for (size_t j = 0; j < 2; ++j) {
            size_t nposition_restraints = file_.read_single_size_as_i32();
            file_.skip(nposition_restraints * 3 * header_.sizeof_real());
        }
    }

    const size_t natoms = file_.read_single_size_as_i32();
    assert(natoms == header_.natoms);

    if (header_.file_version >= TPRVersion::IntermolecularBondeds) {
        bool has_intermolecular_bonds = read_gmx_bool();
        if (has_intermolecular_bonds) {
            InteractionLists interaction_lists =
                read_interaction_lists(file_, header_.file_version);
            add_conectivity(frame, interaction_lists);
        }
    }

    // Skip atom types for old formats
    // see `do_atomtypes`
    if (header_.file_version < TPRVersion::RemoveAtomtypes) {
        size_t ntypes = file_.read_single_size_as_i32();
        if (header_.file_version < TPRVersion::RemoveImplicitSolvation) {
            file_.skip(3 * ntypes * header_.sizeof_real());
        }
        file_.skip(ntypes * sizeof(int32_t));
        if (header_.file_version >= 60 &&
            header_.file_version < TPRVersion::RemoveImplicitSolvation) {
            file_.skip(2 * ntypes * header_.sizeof_real());
        }
    }

    // Skip dihedral correction maps (CMAP)
    // see `do_cmap`
    if (header_.file_version >= 65) {
        size_t ngrids = file_.read_single_size_as_i32();
        size_t grid_spacing = file_.read_single_size_as_i32();
        file_.skip(ngrids * grid_spacing * grid_spacing * 4 * header_.sizeof_real());
    }

    // Skip atom groups
    // see `do_groups`
    {
        // Skip sizes of groups
        for (size_t i = 0; i < NR_GROUP_TYPES; ++i) {
            // see `do_grps`
            size_t group_size = file_.read_single_size_as_i32();
            file_.skip(group_size * sizeof(int32_t));
        }
        size_t ngroup_names = file_.read_single_size_as_i32();
        // Skip symbol table indices
        file_.skip(ngroup_names * sizeof(int32_t));
        // Skip group numbers
        for (size_t i = 0; i < NR_GROUP_TYPES; ++i) {
            size_t ngroup_numbers = file_.read_single_size_as_i32();
            if (header_.body_convention == FileIOXdr) {
                file_.skip(ngroup_numbers * sizeof(uint32_t));
            } else {
                assert(header_.body_convention == InMemory);
                file_.skip(ngroup_numbers * sizeof(uint8_t));
            }
        }
    }

    if (header_.file_version >= TPRVersion::StoreNonBondedInteractionExclusionGroup) {
        int64_t intermolecularExclusionGroupSize = file_.read_single_i64();
        if (intermolecularExclusionGroupSize < 0) {
            throw format_error("invalid intermolecular exclusion group size in TPR file: expected "
                               "a positive value, got {}",
                               intermolecularExclusionGroupSize);
        }
        uint64_t skip_size =
            static_cast<uint64_t>(intermolecularExclusionGroupSize) * sizeof(int32_t);
        file_.skip(skip_size);
    }
}

void TPRFormat::read_coordinates(Frame& frame) {
    if (header_.use_double) {
        // Double
        std::vector<double> dx(header_.natoms * 3);
        if (header_.has_positions) {
            file_.read_f64(dx);
            auto positions = frame.positions();
            assert(dx.size() == 3 * positions.size());
            for (size_t i = 0; i < frame.size(); i++) {
                // Factor 10 because the cell lengths are in nm in the TPR format
                positions[i][0] = dx[i * 3] * 10.0;
                positions[i][1] = dx[i * 3 + 1] * 10.0;
                positions[i][2] = dx[i * 3 + 2] * 10.0;
            }
        }
        if (header_.has_velocities) {
            file_.read_f64(dx);
            frame.add_velocities();
            auto velocities = *frame.velocities();
            assert(dx.size() == 3 * velocities.size());
            for (size_t i = 0; i < frame.size(); i++) {
                // Factor 10 because the lengths are in nm in the TPR format
                velocities[i][0] = dx[i * 3] * 10.0;
                velocities[i][1] = dx[i * 3 + 1] * 10.0;
                velocities[i][2] = dx[i * 3 + 2] * 10.0;
            }
        }
    } else {
        // Float
        std::vector<float> dx(header_.natoms * 3);
        if (header_.has_positions) {
            file_.read_f32(dx);
            auto positions = frame.positions();
            assert(dx.size() == 3 * positions.size());
            for (size_t i = 0; i < frame.size(); i++) {
                // Factor 10 because the cell lengths are in nm in the TPR format
                positions[i][0] = static_cast<double>(dx[i * 3]) * 10.0;
                positions[i][1] = static_cast<double>(dx[i * 3 + 1]) * 10.0;
                positions[i][2] = static_cast<double>(dx[i * 3 + 2]) * 10.0;
            }
        }
        if (header_.has_velocities) {
            file_.read_f32(dx);
            frame.add_velocities();
            auto velocities = *frame.velocities();
            assert(dx.size() == 3 * velocities.size());
            for (size_t i = 0; i < frame.size(); i++) {
                // Factor 10 because the lengths are in nm in the TPR format
                velocities[i][0] = static_cast<double>(dx[i * 3]) * 10.0;
                velocities[i][1] = static_cast<double>(dx[i * 3 + 1]) * 10.0;
                velocities[i][2] = static_cast<double>(dx[i * 3 + 2]) * 10.0;
            }
        }
    }
    if (header_.has_forces) {
        file_.skip(header_.natoms * 3 * header_.sizeof_real());
    }
}

std::vector<std::string> TPRFormat::read_symbol_table() {
    size_t nsymbols = file_.read_single_size_as_i32();
    std::vector<std::string> symbol_table;
    symbol_table.reserve(nsymbols);
    for (size_t i = 0; i < nsymbols; ++i) {
        symbol_table.emplace_back(read_gmx_string());
    }
    return symbol_table;
}

const std::string& TPRFormat::read_symbol_table_entry(const std::vector<std::string>& table) {
    size_t idx = file_.read_single_size_as_i32();
    return table.at(idx);
}

std::string TPRFormat::read_gmx_string() {
    if (header_.body_convention == FileIOXdr) {
        return file_.read_gmx_string();
    } else {
        assert(header_.body_convention == InMemory);
        // Read a non-XDR-compliant *long* GROMACS string
        // A long GROMACS string stores the length of the string as uint64 before
        // the string contents. The contents are *not* padded to make the size a
        // multiple of four bytes.

        // lenght without null terminator
        const uint64_t len = file_.read_single_u64();
        // next comes the contents of the string without padding
        std::vector<char> buf;
        buf.resize(static_cast<size_t>(len));
        file_.read_char(buf);
        return std::string(buf.begin(), buf.end());
    }
}

uint8_t TPRFormat::read_gmx_uchar() {
    if (header_.body_convention == FileIOXdr) {
        return static_cast<uint8_t>(file_.read_single_u32());
    } else {
        assert(header_.body_convention == InMemory);
        return file_.read_single_u8();
    }
}

bool TPRFormat::read_gmx_bool() {
    if (header_.body_convention == FileIOXdr) {
        return file_.read_single_i32() != 0;
    } else {
        assert(header_.body_convention == InMemory);
        return file_.read_single_u8() != 0;
    }
}
