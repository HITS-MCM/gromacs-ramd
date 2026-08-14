/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#pragma once

#include <random>

#include "gromacs/random/seed.h"
#include "gromacs/random/threefry.h"
#include "gromacs/random/uniformrealdistribution.h"
#include "gromacs/utility/vectypes.h"

namespace gmx
{

/*! \internal \brief
 * Common interface for generating random points on the unit sphere, so that
 * RAMD can switch between RNG implementations at run time (see
 * RandomSphericalDirectionGenerator and LegacyRandomSphericalDirectionGenerator).
 */
class IRandomSphericalDirectionGenerator
{
public:
    virtual ~IRandomSphericalDirectionGenerator() = default;

    virtual DVec operator()() = 0;
};

class RandomSphericalDirectionGenerator final : public IRandomSphericalDirectionGenerator
{
public:
    // ThreeFry2x64, unlike std::default_random_engine, is fully specified by GROMACS
    // rather than left to the standard library implementation, so a given seed
    // produces the same sequence of directions regardless of compiler/platform.
    RandomSphericalDirectionGenerator(int64_t seed) : engine_(seed, RandomDomain::Other) {}

    DVec operator()() override
    {
        // azimuth angle
        real theta = 2 * M_PI * dist_(engine_);

        // polar angle
        real psi = std::acos(1.0 - 2 * dist_(engine_));

        DVec direction;
        direction[0] = std::cos(theta) * std::sin(psi);
        direction[1] = std::sin(theta) * std::sin(psi);
        direction[2] = std::cos(psi);

        return direction;
    }

private:
    /// Random number generator
    ThreeFry2x64<> engine_;

    /// Random number distribution
    UniformRealDistribution<real> dist_;
};

/*! \internal \brief
 * Reproduces the RNG used by RandomSphericalDirectionGenerator prior to the switch
 * to GROMACS' portable ThreeFry2x64 engine. std::default_random_engine's sequence is
 * left to the standard library implementation, so it differs between compilers and
 * platforms; this is kept only so that a ramd-seed from before that switch can still
 * reproduce the exact same trajectory (via the ramd-legacy-rng mdp option).
 */
class LegacyRandomSphericalDirectionGenerator final : public IRandomSphericalDirectionGenerator
{
public:
    LegacyRandomSphericalDirectionGenerator(int64_t seed) : engine_(seed), dist_(0.0, 1.0) {}

    DVec operator()() override
    {
        // azimuth angle
        real theta = 2 * M_PI * dist_(engine_);

        // polar angle
        real psi = std::acos(1.0 - 2 * dist_(engine_));

        DVec direction;
        direction[0] = std::cos(theta) * std::sin(psi);
        direction[1] = std::sin(theta) * std::sin(psi);
        direction[2] = std::cos(psi);

        return direction;
    }

private:
    /// Random number generator
    std::default_random_engine engine_;

    /// Random number distribution
    std::uniform_real_distribution<> dist_;
};

} // namespace gmx
