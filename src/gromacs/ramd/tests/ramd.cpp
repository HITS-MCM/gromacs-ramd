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

TEST(RAMDTest, construction)
{
    RAMDParams params;
    RAMD ramd(params);
}

} // namespace test
} // namespace gmx
