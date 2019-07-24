/*
 * RAMDTest.cpp
 *
 *  Created on: Jul 24, 2019
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <gtest/gtest.h>
#include <iostream>
#include <string>
#include <vector>

#include "gromacs/options/filenameoption.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/real.h"
#include "gromacs/gmxpreprocess/grompp.h"
#include "programs/mdrun/mdrun_main.h"
#include "testutils/cmdlinetest.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{
namespace
{

struct TestDataStructure
{
    TestDataStructure(std::string const& testDirectory)
      : testDirectory(testDirectory)
    {}

    std::string testDirectory;
};

//! Test fixture for RAMD
class RAMDTest : public ::testing::WithParamInterface<TestDataStructure>,
                public CommandLineTestBase
{};

//! Test body for RAMD
TEST_P(RAMDTest, Basic)
{
    std::cout << GetParam().testDirectory << std::endl;

    std::string cwd = gmx::Path::getWorkingDirectory();
    std::string dataPath = std::string(fileManager().getInputDataDirectory()) + "/data";
    std::string testPath = fileManager().getTemporaryFilePath("/" + GetParam().testDirectory);

    std::string cmd = "mkdir -p " + testPath;
    ASSERT_FALSE(system(cmd.c_str()));

    cmd = "cp -r " + dataPath + "/" + GetParam().testDirectory + "/* " + testPath;
    ASSERT_FALSE(system(cmd.c_str()));

    gmx_chdir(testPath.c_str());

    ::gmx::test::CommandLine call_grompp;
    call_grompp.append("gmx grompp");
    call_grompp.addOption("-f", "ramd.mdp");

    std::cout << "grompp command: " << call_grompp.toString() << std::endl;

    ::gmx::test::CommandLine call_mdrun;
    call_mdrun.append("gmx mdrun");
    call_mdrun.addOption("-s", "topol.tpr");

    std::cout << "mdrun command: " << call_mdrun.toString() << std::endl;
}

INSTANTIATE_TEST_CASE_P(DISABLED_AllRAMDTests, RAMDTest, ::testing::Values(
    TestDataStructure("1WDHI"),
    TestDataStructure("membrane")
));

} // namespace
} // namespace test
} // namespace gmx
