/*
 * RandomSphericalDirectionGenerator.h
 *
 *  Created on: Jun 25, 2019
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#pragma once

#include <iostream>
#include <random>

#include "gromacs/math/vectypes.h"

namespace gmx {

class RandomSphericalDirectionGenerator
{
public:

	RandomSphericalDirectionGenerator(int64_t seed, bool use_old_angle_dist = false)
     : engine(seed),
       dist(0.0, 1.0),
	   use_old_angle_dist(use_old_angle_dist)
    {
		if (use_old_angle_dist) {
			std::cout << "==== RAMD ==== Warning: Old angle distribution is used." << std::endl;
		}
    }

	DVec operator () ()
	{
	    // azimuth angle
	    real theta = 2 * M_PI * dist(engine);

	    // polar angle
	    real psi;
	    if (use_old_angle_dist) psi = M_PI * dist(engine);
	    else psi = std::acos(1.0 - 2 * dist(engine));

	    DVec direction;
	    direction[0] = std::cos(theta) * std::sin(psi);
	    direction[1] = std::sin(theta) * std::sin(psi);
	    direction[2] = std::cos(psi);

	    return direction;
	}

private:

	/// Random number generator
	std::default_random_engine engine;

	/// Random number distribution
    std::uniform_real_distribution<> dist;

    /// For backward compa
    bool use_old_angle_dist;

};

} // namespace gmx
