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

	RAMD(RAMDParams const& params, pull_t *pull);

	real add_force(int64_t step, t_mdatoms const& mdatoms,
		gmx::ForceWithVirial *forceWithVirial);

private:

	/// Returns a new random direction
	void set_random_direction() const;

	/// Initialization number for pseudo random number generator
	const RAMDParams params;

	/// Pointer to pull work array
	pull_t *pull;

	/// Random number generator
	mutable std::default_random_engine engine;

	/// Random number distribution
    mutable std::uniform_real_distribution<> dist;

    /// COM of receptor of last RAMD evaluation step
    DVec com_rec_prev;

    /// COM of ligand of last RAMD evaluation step
    DVec com_lig_prev;

};

} // namespace gmx
