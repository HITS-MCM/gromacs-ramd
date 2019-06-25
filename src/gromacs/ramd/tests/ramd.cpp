/**
 * \internal \file
 *
 * \brief
 * Implements test of RAMD
 *
 * \author Bernd Doser <bernd.doser@h-its.org>
 *
 */

#include <gtest/gtest.h>

#include "gromacs/mdtypes/state.h"
#include "gromacs/pulling/pull_internal.h"
#include "gromacs/ramd/ramd.h"

namespace gmx {
namespace test {

TEST(RAMDTest, exception_group_number)
{
	RAMDParams params;
	pull_t pull;

	EXPECT_THROW(RAMD(params, &pull), std::exception);
}

TEST(RAMDTest, construction)
{
	RAMDParams params;
	pull_t pull;
//	pull.group.resize(3);
//
	EXPECT_THROW(RAMD(params, &pull), std::exception);
}

} // namespace test
} // namespace gmx
