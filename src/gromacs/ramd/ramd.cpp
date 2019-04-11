/*
 * ramd.cpp
 *
 *  Created on: Apr 10, 2019
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include "ramd.h"
#include "gromacs/mdtypes/pull_params.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/pulling/pull_internal.h"

real RAMD::add_force(struct pull_t *pull,
	                 const t_mdatoms &mdatoms,
                     gmx::ForceWithVirial *forceWithVirial) const
{
	pull->coord[0].spatialData.vec[0] = 0.0;
	pull->coord[0].spatialData.vec[1] = 0.0;
	pull->coord[0].spatialData.vec[2] = 1.0;

	double force = 10.0;

    apply_external_pull_coord_force(pull, 0, force, &mdatoms, forceWithVirial);

    return 0.0;
}
