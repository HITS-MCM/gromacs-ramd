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
}

// Pseudo-Code:
//{
//	set_random_direction();
//
//	if (step % params.ramd_eval_freq) {
//
//	}
//}

real RAMD::add_force(t_mdatoms const& mdatoms,
                     gmx::ForceWithVirial *forceWithVirial) const
{
	set_random_direction(pull->coord[0].spatialData.vec);

	double force = 10.0;

    apply_external_pull_coord_force(pull, 0, force, &mdatoms, forceWithVirial);

    return 0.0;
}

void RAMD::set_random_direction(dvec& vec) const
{
//	set theta [expr "2*$pi*rand()"]
//	set psi [expr "$pi*rand()"]
//	set sinpsi [expr "sin($psi)"]
//	set rx [expr "cos($theta)*$sinpsi"]
//	set ry [expr "sin($theta)*$sinpsi"]
//	set rz [expr "cos($psi)"]

    auto r = dist(engine);
	vec[0] = std::cos(r);
	vec[1] = std::sin(r);
	vec[2] = std::cos(r);
}

} // namespace gmx
