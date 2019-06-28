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
#include <vector>

#include "gromacs/math/vectypes_io.h"
#include "gromacs/ramd/RandomSphericalDirectionGenerator.h"

namespace gmx {
namespace test {

std::vector<int> create_histogram(int number_of_buckets, std::vector<double> const& data)
{
	const double bucket_size = 2 * M_PI / number_of_buckets;
	std::vector<int> histogram(number_of_buckets, 0);

	for (auto const& d : data) ++histogram[floor((d + M_PI) / bucket_size)];

	return histogram;
}

void check_distribution(int x, int y, std::vector<DVec> const& random_directions)
{
	std::vector<double> thetas;
	for (auto const& d : random_directions)
	{
		thetas.push_back(std::atan2(d[y], d[x]));
	}

	auto histogram = create_histogram(32, thetas);

	for (auto const& h : histogram)
	{
		std::cout << h << std::endl;
	}
}

TEST(RandomSphericalDirectionGeneratorTest, angle_distribution)
{
    RandomSphericalDirectionGenerator random_spherical_direction_generator(1234);

    std::vector<DVec> random_directions;
    for (int i = 0; i < 1000000; ++i)
    {
    	random_directions.push_back(random_spherical_direction_generator());
    }

    for (auto const& d : random_directions)
    {
        EXPECT_NEAR(1.0, d.norm2(), 1e-6);
    }

	std::cout << "xy" << std::endl;
    check_distribution(0, 1, random_directions);

	std::cout << "xz" << std::endl;
	check_distribution(0, 2, random_directions);

	std::cout << "zy" << std::endl;
	check_distribution(2, 1, random_directions);
}

} // namespace test
} // namespace gmx
