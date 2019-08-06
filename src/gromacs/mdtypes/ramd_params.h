#pragma once

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

namespace gmx {

/// Parameters for RAMD group
struct RAMDGroup
{
    int   nat;     ///< Number of atoms in the RAMD group
    int  *ind;     ///< The global atoms numbers
};

/// Parameters for Random Acceleration Molecular Dynamics (RAMD)
struct RAMDParams
{
    int64_t        seed;             ///< Initialization number for pseudo random number generator
    real           force;            ///< Force to be applied in kcal/mol/Angstrom
    int            eval_freq;        ///< Number of MD steps in one RAMD step
    real           r_min_dist;       ///< Specifies the minimum distance in Angstrom
                                     ///  to be traveled by the ligand in one RAMD step
    int            force_out_freq;   ///< Every 'force_out_freq' steps detailed output of forces will be written
    real           max_dist;         ///< Specifies the distance in Angstrom between the COMs of the ligand
                                     ///  and the protein when the simulation is stopped
    RAMDGroup      protein;          ///< The name of the protein group, is looked up in the index file or
                                     ///  in the default groups to obtain the atoms involved
    RAMDGroup      ligand;           ///< The name of the protein group, is looked up in the index file or
                                     ///  in the default groups to obtain the atoms involved
    gmx_bool       old_angle_dist;   ///< Use old angle distribution
};

} // namespace gmx
