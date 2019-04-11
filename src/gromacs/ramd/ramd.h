/*
 * ramd.h
 *
 *  Created on: Apr 10, 2019
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#pragma once

#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/forceoutput.h"

class RAMD
{
public:

	RAMD(int seed)
     : seed(seed)
    {}

	real add_force(struct pull_t *pull,
			       const t_mdatoms &mdatoms,
                   gmx::ForceWithVirial *forceWithVirial) const;

private:

	/// Initialization number for pseudo random number generator
	int seed;

};
