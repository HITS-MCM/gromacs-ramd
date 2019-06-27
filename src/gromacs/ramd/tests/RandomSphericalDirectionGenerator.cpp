/**
 * \internal \file
 *
 * \brief
 * Implements test of RandomSphericalDirectionGenerator
 *
 * \author Bernd Doser <bernd.doser@h-its.org>
 *
 */

#include <gtest/gtest.h>

#include "gromacs/math/vectypes_io.h"
#include "gromacs/ramd/RandomSphericalDirectionGenerator.h"

namespace gmx {
namespace test {

TEST(RandomSphericalDirectionGeneratorTest, angle_distribution)
{
    RandomSphericalDirectionGenerator random_spherical_direction_generator(1234);
    DVec d = random_spherical_direction_generator();
    std::cout << "d = " << d << std::endl;
}

} // namespace test
} // namespace gmx
