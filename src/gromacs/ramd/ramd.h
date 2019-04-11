/*
 * ramd.h
 *
 *  Created on: Apr 10, 2019
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#pragma once

#include <random>

#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/ramd_params.h"
#include "gromacs/mdtypes/forceoutput.h"

struct pull_t;

namespace gmx {

class RAMD
{
public:

	RAMD(RAMDParams const* params);

	real add_force(pull_t *pull,
			       const t_mdatoms &mdatoms,
                   gmx::ForceWithVirial *forceWithVirial) const;

private:

	/// Returns a new random direction
	void set_random_direction(dvec& vec) const;

	/// Initialization number for pseudo random number generator
	RAMDParams const* params;

	/// Random number generator
	mutable std::default_random_engine engine;

	/// Random number distribution
    mutable std::uniform_real_distribution<> dist;

};

} // namespace gmx
