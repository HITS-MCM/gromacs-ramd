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
#include "gmxpre.h"

#include "ramd.h"

#include <cassert>

#include "gromacs/commandline/filenm.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/mdlib/sighandler.h"
#include "gromacs/mdtypes/pull_params.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/pulling/pull_internal.h"
#include "gromacs/utility/fatalerror.h"

namespace gmx
{

RAMD::RAMD(const RAMDParams&           params,
           pull_t*                     pull,
           int64_t*                    pstep,
           const gmx::StartingBehavior startingBehavior,
           const t_commrec*            cr,
           int                         nfile,
           const t_filenm              fnm[],
           const gmx_output_env_t*     oenv)
  : params(params),
    pull(pull),
    pstep(pstep),
    random_spherical_direction_generator(params.seed, params.old_angle_dist),
    direction(params.ngroup),
    com_rec_prev(params.ngroup),
    com_lig_prev(params.ngroup),
    out(nullptr),
    cr(cr)
{
    for (int g = 0; g < params.ngroup; ++g)
    {
        register_external_pull_potential(pull, g * 3    , "RAMD");
        register_external_pull_potential(pull, g * 3 + 1, "RAMD");
        register_external_pull_potential(pull, g * 3 + 2, "RAMD");
        direction[g] = random_spherical_direction_generator();
    }

    if (MASTER(cr) and opt2bSet("-ramd", nfile, fnm))
    {
        auto filename = std::string(opt2fn("-ramd", nfile, fnm));

        if (startingBehavior == gmx::StartingBehavior::RestartWithAppending)
        {
            out = gmx_fio_fopen(filename.c_str(), "a+");
        }
        else
        {
            out = gmx_fio_fopen(filename.c_str(), "w+");
            xvgr_header(out, "RAMD Ligand-receptor COM distance", "Time (ps)", "Distance (nm)", exvggtXNY, oenv);
        }
        fflush(out);
    }
}

void RAMD::calculateForces(const ForceProviderInput& forceProviderInput,
                           ForceProviderOutput*      forceProviderOutput)
{
    if (*pstep == 0)
    {
        // Store COM positions for first evaluation
        for (int g = 0; g < params.ngroup; ++g)
        {
            com_rec_prev[g] = pull->group[g * 2].x;
            com_lig_prev[g] = pull->group[g * 2 + 1].x;

            if (MASTER(cr) and out)
            {
                auto dist = std::sqrt((com_lig_prev[g] - com_rec_prev[g]).norm2());
                fprintf(out, "%.4f\t%g\n", forceProviderInput.t_, dist);
                fflush(out);
            }
        }
    }
    else if (!(*pstep % params.eval_freq))
    {
        if (MASTER(cr) and debug)
        {
            fprintf(debug, "==== RAMD ==== evaluation %ld\n", *pstep);
        }
        for (int g = 0; g < params.ngroup; ++g)
        {
            DVec com_rec_curr = pull->group[g * 2].x;
            DVec com_lig_curr = pull->group[g * 2 + 1].x;
            auto curr_dist = std::sqrt((com_lig_curr - com_rec_curr).norm2());

            if (MASTER(cr) and debug)
            {
                fprintf(debug, "==== RAMD ==== group %d\n", g);
                fprintf(debug, "==== RAMD ==== COM ligand position at [%g, %g, %g]\n",
                        com_lig_curr[0], com_lig_curr[1], com_lig_curr[2]);
                fprintf(debug, "==== RAMD ==== COM receptor position at [%g, %g, %g]\n",
                        com_rec_curr[0], com_rec_curr[1], com_rec_curr[2]);
                fprintf(debug,
                        "==== RAMD ==== Distance between COM of receptor and COM of ligand is %g\n",
                        curr_dist);
            }

            if (MASTER(cr) and out)
            {
                fprintf(out, "%.4f\t%g\n", forceProviderInput.t_, curr_dist);
                fflush(out);
            }

            if (MASTER(cr) and curr_dist >= params.group[0].max_dist)
            {
                if (debug)
                {
                    fprintf(debug,
                            "==== RAMD ==== Maximal distance between ligand and receptor COM is "
                            "reached.\n");
                }
                fprintf(stdout, "==== RAMD ==== GROMACS will be stopped after %ld steps.\n", *pstep);
                gmx_set_stop_condition(gmx_stop_cond_next);
            }

            // walk_dist = vector length of the vector substraction
            // (com_lig_curr - com_rec_curr) - (com_lig_prev - com_rec_prev)
            // during last RAMD evaluation step
            auto walk_dist =
                    std::sqrt(((com_lig_curr - com_rec_curr) - (com_lig_prev[g] - com_rec_prev[g])).norm2());

            if (MASTER(cr) and debug)
            {
                fprintf(debug, "==== RAMD ==== Previous COM ligand position at [%g, %g, %g]\n",
                        com_lig_prev[g][0], com_lig_prev[g][1], com_lig_prev[g][2]);
                fprintf(debug, "==== RAMD ==== Previous COM receptor position at [%g, %g, %g]\n",
                        com_rec_prev[g][0], com_rec_prev[g][1], com_rec_prev[g][2]);
                fprintf(debug,
                        "==== RAMD ==== Change in receptor-ligand distance since last RAMD evaluation "
                        "is %g\n",
                        walk_dist);
            }

            if (walk_dist < params.group[0].r_min_dist)
            {
                direction[g] = random_spherical_direction_generator();
                if (MASTER(cr) and debug)
                {
                    fprintf(debug, "==== RAMD ==== New random direction is [%g, %g, %g]\n",
                            direction[g][0], direction[g][1], direction[g][2]);
                }
            }

            com_lig_prev[g] = com_lig_curr;
            com_rec_prev[g] = com_rec_curr;
        }
    }

    t_pbc pbc;
    set_pbc(&pbc, pull->ePBC, forceProviderInput.box_);

    for (int g = 0; g < params.ngroup; ++g)
    {
        for (int i = 0; i < 3; ++i)
        {
            get_pull_coord_value(pull, g * 2 + i, &pbc);
            apply_external_pull_coord_force(pull, g * 2 + i, direction[g][i] * params.group[g].force,
                &forceProviderInput.mdatoms_, &forceProviderOutput->forceWithVirial_);
        }
    }
}

} // namespace gmx
