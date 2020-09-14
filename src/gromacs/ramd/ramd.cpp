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

#include "gromacs/mdtypes/pull_params.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/pulling/pull_internal.h"
#include "gromacs/utility/fatalerror.h"

namespace gmx
{

RAMD::RAMD(RAMDParams const& params) :
    params(params),
    random_spherical_direction_generator(params.seed, params.old_angle_dist),
    direction(random_spherical_direction_generator())
{
}

real RAMD::add_force(int64_t               step,
                     t_mdatoms const&      mdatoms,
                     gmx::ForceWithVirial* forceWithVirial,
                     pull_t*               pull,
                     const t_commrec*      cr)
{
    assert(pull->group.size() == 3);

    if (MASTER(cr))
    {
        if (step == 0)
        {
            // Store COM positions for first evaluation
            com_rec_prev = pull->group[1].x;
            com_lig_prev = pull->group[2].x;
        }
        else if (!(step % params.eval_freq))
        {
            fprintf(debug, "==== RAMD ==== evaluation %ld\n", step);

            DVec com_rec_curr = pull->group[1].x;
            DVec com_lig_curr = pull->group[2].x;

            fprintf(debug, "==== RAMD ==== COM ligand position at [%g, %g, %g]\n", com_lig_curr[0],
                    com_lig_curr[1], com_lig_curr[2]);
            fprintf(debug, "==== RAMD ==== COM receptor position at [%g, %g, %g]\n",
                    com_rec_curr[0], com_rec_curr[1], com_rec_curr[2]);

            auto curr_dist = std::sqrt((com_lig_curr - com_rec_curr).norm2());

            fprintf(debug,
                    "==== RAMD ==== Distance between COM of receptor and COM of ligand is %g\n",
                    curr_dist);

            if (curr_dist >= params.max_dist)
            {
                fprintf(debug,
                        "==== RAMD ==== Maximal distance between ligand and receptor COM is "
                        "reached.\n");
                fprintf(stdout, "==== RAMD ==== GROMACS will be stopped after %ld steps.\n", step);
                std::abort();
            }

            fprintf(debug, "==== RAMD ==== Previous COM ligand position at [%g, %g, %g]\n",
                    com_lig_prev[0], com_lig_prev[1], com_lig_prev[2]);
            fprintf(debug, "==== RAMD ==== Previous COM receptor position at [%g, %g, %g]\n",
                    com_rec_prev[0], com_rec_prev[1], com_rec_prev[2]);

            // walk_dist = vector length of the vector substraction
            // (com_lig_curr - com_rec_curr) - (com_lig_prev - com_rec_prev)
            // during last RAMD evaluation step
            auto walk_dist =
                    std::sqrt(((com_lig_curr - com_rec_curr) - (com_lig_prev - com_rec_prev)).norm2());

            fprintf(debug,
                    "==== RAMD ==== Change in receptor-ligand distance since last RAMD evaluation "
                    "is %g\n",
                    walk_dist);

            if (walk_dist < params.r_min_dist)
            {
                direction = random_spherical_direction_generator();
                fprintf(debug, "==== RAMD ==== New random direction is [%g, %g, %g]\n",
                        direction[0], direction[1], direction[2]);
            }

            com_lig_prev = com_lig_curr;
            com_rec_prev = com_rec_curr;
        }
    }

    real potential = 0.0;
    apply_external_pull_coord_force(pull, 0, direction[0] * params.force, &mdatoms, forceWithVirial);
    apply_external_pull_coord_force(pull, 1, direction[1] * params.force, &mdatoms, forceWithVirial);
    apply_external_pull_coord_force(pull, 2, direction[2] * params.force, &mdatoms, forceWithVirial);
    return potential;
}

} // namespace gmx
