#ifndef GMX_MDTYPES_RAMD_PARAMS_H
#define GMX_MDTYPES_RAMD_PARAMS_H

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

namespace gmx {

//! Parameters for RAMD group
struct RAMDGroup
{
    int   nat;     ///< Number of atoms in the RAMD group
    int  *ind;     ///< The global atoms numbers
};

//! Parameters for Random Acceleration Molecular Dynamics (RAMD)
struct RAMDParams
{
    gmx_int64_t    seed;           ///< Random seed
    real           force;          ///< Force to be applied in kcal/mol/Angstrom
    int            steps;          ///< Number of MD steps in one RAMD step
    real           rMinRamd;       ///< Specifies the minimum distance in Angstrom
                                   ///  to be traveled by the ligand in one RAMD step
    int            forceOutFreq;   ///< Every 'forceOutFreq' steps detailed output of forces will be written
    real           maxDist;        ///< Specifies the distance in Angstrom between the COMs of the ligand
                                   ///  and the protein when the simulation is stopped
    RAMDGroup      protein;        ///< The name of the protein group, is looked up in the index file or
                                   ///  in the default groups to obtain the atoms involved
    RAMDGroup      ligand;         ///< The name of the protein group, is looked up in the index file or
                                   ///  in the default groups to obtain the atoms involved
};

} // namespace gmx

#endif /* GMX_MDTYPES_RAMD_PARAMS_H */
