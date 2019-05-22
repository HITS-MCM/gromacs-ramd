/*
 * ramd.cpp
 *
 *  Created on: Apr 10, 2019
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <iostream>

#include "ramd.h"
#include "gromacs/mdtypes/pull_params.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/pulling/pull_internal.h"

namespace gmx {

RAMD::RAMD(RAMDParams const& params, pull_t *pull)
 : params(params),
   pull(pull),
   engine(params.seed),
   dist(0.0, 1.0),
   // Store COM positions for first evaluation
   com_rec_prev(pull->group[1].x),
   com_lig_prev(pull->group[2].x)
{
	set_random_direction();
}

real RAMD::add_force(int64_t step, t_mdatoms const& mdatoms,
    gmx::ForceWithVirial *forceWithVirial)
{
    if (step != 0 && !(step % params.eval_freq))
	{
	    std::cout << "==== RAMD ==== evaluation " << step << std::endl;

	    DVec com_rec_curr = pull->group[1].x;
	    DVec com_lig_curr = pull->group[2].x;

	    std::cout << "==== RAMD ==== COM ligand position at [" << com_lig_curr[0] << " ," << com_lig_curr[1] << " ," << com_lig_curr[2] << "]" << std::endl;
	    std::cout << "==== RAMD ==== COM receptor position at [" << com_rec_curr[0] << " ," << com_rec_curr[1] << " ," << com_rec_curr[2] << "]" << std::endl;

        auto curr_dist = std::sqrt((com_lig_curr - com_rec_curr).norm2());

        std::cout << "==== RAMD ==== Distance between COM of receptor and COM of ligand is " << curr_dist << std::endl;

        if (curr_dist >= params.max_dist) {
        	std::cout << "==== RAMD ==== Maximal distance between ligand and receptor COM is reached.\n"
        			  << "==== RAMD ==== GROMACS will be stopped." << std::endl;
            std::abort();
        }

	    std::cout << "==== RAMD ==== Previous COM ligand position at " << com_lig_prev[0] << " ," << com_lig_prev[1] << " ," << com_lig_prev[2] << "]" << std::endl;
	    std::cout << "==== RAMD ==== Previous COM receptor position at " << com_rec_prev[0] << " ," << com_rec_prev[1] << " ," << com_rec_prev[2] << "]" << std::endl;

        // walk_dist =	vector length of the vector substraction
        // (com_lig_curr - com_rec_curr) - (com_lig_prev - com_rec_prev)
        // during last RAMD evaluation step
        auto walk_dist = std::sqrt(((com_lig_curr - com_rec_curr) - (com_lig_prev - com_rec_prev)).norm2());

        std::cout << "==== RAMD ==== Change in receptor-ligand distance since last RAMD evaluation is " << walk_dist << std::endl;

        if (walk_dist < params.r_min_dist) {
            set_random_direction();
        }

        com_lig_prev = com_lig_curr;
        com_rec_prev = com_rec_curr;
	}

	real potential = 0.0;
    apply_external_pull_coord_force(pull, 0, params.force, &mdatoms, forceWithVirial);
    return potential;
}

void RAMD::set_random_direction() const
{
    auto theta = 2 * M_PI * dist(engine);
    auto psi   =     M_PI * dist(engine);

    auto& vec = pull->coord[0].spatialData.vec;
	vec[0] = std::cos(theta) * std::sin(psi);
	vec[1] = std::sin(theta) * std::sin(psi);
	vec[2] = std::cos(psi);

    std::cout << "==== RAMD ==== New random direction is [" << vec[0] << " ," << vec[1] << " ," << vec[2] << "]" << std::endl;
}

} // namespace gmx
