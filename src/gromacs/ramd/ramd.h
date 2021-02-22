/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017 by the GROMACS development team.
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#pragma once

#include <memory>

#include "gromacs/fileio/oenv.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/ramd_params.h"
#include "gromacs/ramd/randomsphericaldirectiongenerator.h"

struct pull_t;

namespace gmx
{

class RAMD final : public IForceProvider
{
public:

    /*! \brief RAMD ForceProvider.
    *
    * \param[in]     ramdParams              RAMD parameters.
    * \param[in,out] pull                    Pointer to pull object.
    * \param[in]     startingBehavior        Describes whether this is a restart appending to output
    *                                        files.
    * \param[in]     cr                      Struct for communication, can be nullptr.
    * \param[in]     nfile                   Number of files.
    * \param[in]     fnm                     Filename struct.
    * \param[in]     oenv                    The output environment information.
    */
    RAMD(const RAMDParams&           params,
         pull_t*                     pull,
         int64_t*                    pstep,
         const gmx::StartingBehavior startingBehavior,
         const t_commrec*            cr,
         int                         nfile,
         const t_filenm              fnm[],
         const gmx_output_env_t*     oenv);

    //! \copydoc IForceProvider::calculateForces()
    void calculateForces(const ForceProviderInput& forceProviderInput,
                         ForceProviderOutput*      forceProviderOutput) override;

private:

    /// RAMD parameters
    const RAMDParams& params;

    /// PULL work structure
    pull_t* pull;

    /// Pointer to the current MD step
    int64_t* pstep;

    /// Random pull direction
    RandomSphericalDirectionGenerator random_spherical_direction_generator;

    /// Current pull direction
    std::vector<DVec> direction;

    /// COM of receptor of last RAMD evaluation step
    std::vector<DVec> com_rec_prev;

    /// COM of ligand of last RAMD evaluation step
    std::vector<DVec> com_lig_prev;

    /// Output file for COM distances
    FILE* out;

    /// MPI communicator
    const t_commrec* cr;

    /// Has the ligand left his binding site?
    std::vector<int> ligand_exited;

};

} // namespace gmx
