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
   dist(0.0, 1.0)
{
    std::cout << "==== RAMD seed = " << params.seed << std::endl;
	set_random_direction();
}

real RAMD::add_force(int64_t step, t_mdatoms const& mdatoms,
    gmx::ForceWithVirial *forceWithVirial) const
{
    std::cout << "==== RAMD step " << step << std::endl;

	if (step % params.eval_freq)
	{

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
}

} // namespace gmx
