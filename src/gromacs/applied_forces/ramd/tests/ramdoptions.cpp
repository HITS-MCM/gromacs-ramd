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
/*! \internal \file
 * \brief
 * Tests for RAMD module options.
 *
 * \author Bernd Doser <bernd.doser@h-its.org>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "gromacs/applied_forces/ramd/ramdoptions.h"

#include <cstdint>

#include <filesystem>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/applied_forces/ramd/ramd.h"
#include "gromacs/mdtypes/imdpoptionprovider_test_helper.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreemdpwriter.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "testutils/testmatchers.h"

namespace gmx
{
namespace
{

//! Content of a minimal, valid RAMD groups file
const char* const c_groupsFileContent =
        "ramd-group {\n"
        "    receptor Protein\n"
        "    ligand Ligand\n"
        "    force 600.0\n"
        "    max-dist 4.0\n"
        "    r-min-dist 0.0025\n"
        "}\n";

//! Build IndexGroupsAndNames containing the groups referenced by c_groupsFileContent
IndexGroupsAndNames ramdIndexGroupsAndNames()
{
    std::vector<IndexGroup> indexGroups;
    indexGroups.push_back({ "Protein", { 0 } });
    indexGroups.push_back({ "Ligand", { 1 } });
    return IndexGroupsAndNames(indexGroups);
}

class RAMDOptionsTest : public ::testing::Test
{
public:
    static KeyValueTreeObject ramdBuildDefaultMdpValues()
    {
        // Prepare MDP inputs
        KeyValueTreeBuilder mdpValueBuilder;
        mdpValueBuilder.rootObject().addValue(std::string(RAMDModuleInfo::sc_name) + "-active",
                                              std::string("true"));
        return mdpValueBuilder.build();
    }

    static KeyValueTreeObject ramdBuildMdpValues()
    {
        // Prepare MDP inputs
        KeyValueTreeBuilder mdpValueBuilder;
        mdpValueBuilder.rootObject().addValue(std::string(RAMDModuleInfo::sc_name) + "-active",
                                              std::string("true"));
        mdpValueBuilder.rootObject().addValue(std::string(RAMDModuleInfo::sc_name) + "-seed",
                                              std::string("42"));
        return mdpValueBuilder.build();
    }
};

TEST_F(RAMDOptionsTest, DefaultParameters)
{
    RAMDOptions ramdOptions;
    const auto  defaultParameters = ramdOptions.parameters();
    EXPECT_FALSE(defaultParameters.active_);
    EXPECT_EQ(1234, defaultParameters.seed_);
    EXPECT_EQ(0, defaultParameters.ngroups_);
}

TEST_F(RAMDOptionsTest, OptionSetsActive)
{
    RAMDOptions ramdOptions;
    test::fillOptionsFromMdpValues(ramdBuildMdpValues(), &ramdOptions);

    EXPECT_TRUE(ramdOptions.active());
    EXPECT_TRUE(ramdOptions.parameters().active_);
    EXPECT_EQ(42, ramdOptions.parameters().seed_);

    // Write parameters to the KVT
    KeyValueTreeBuilder builder;
    MDLogger            logger;
    ramdOptions.setLogger(logger);
    ramdOptions.writeInternalParametersToKvt(builder.rootObject());
    const auto inputTree = builder.build();

    // Retrieve parameters from the KVT
    ramdOptions.readInternalParametersFromKvt(inputTree);

    EXPECT_TRUE(ramdOptions.active());
    EXPECT_TRUE(ramdOptions.parameters().active_);
    EXPECT_EQ(42, ramdOptions.parameters().seed_);
}

TEST_F(RAMDOptionsTest, GroupsFileRelativeToMdpDirectoryIsResolved)
{
    test::TestFileManager       fileManager;
    const std::filesystem::path groupsFilePath =
            fileManager.getTemporaryFilePath("ramd_groups.dat");
    TextWriter::writeFileFromString(groupsFilePath, c_groupsFileContent);

    // Set the groups-file option to just the filename, as if it had been
    // written relative to the .mdp file rather than the current working directory
    KeyValueTreeBuilder mdpValueBuilder;
    mdpValueBuilder.rootObject().addValue(std::string(RAMDModuleInfo::sc_name) + "-active",
                                          std::string("true"));
    mdpValueBuilder.rootObject().addValue(std::string(RAMDModuleInfo::sc_name) + "-groups-file",
                                          groupsFilePath.filename().string());
    RAMDOptions ramdOptionsWithGroupsFile;
    test::fillOptionsFromMdpValues(mdpValueBuilder.build(), &ramdOptionsWithGroupsFile);
    ramdOptionsWithGroupsFile.setMdpFileDirectory(groupsFilePath.parent_path());

    EXPECT_NO_THROW(ramdOptionsWithGroupsFile.setInputGroupIndices(ramdIndexGroupsAndNames()));
    EXPECT_EQ(1, ramdOptionsWithGroupsFile.parameters().ngroups_);
}

TEST_F(RAMDOptionsTest, GroupsFileAbsolutePathIsUnaffectedByMdpDirectory)
{
    test::TestFileManager       fileManager;
    const std::filesystem::path groupsFilePath =
            fileManager.getTemporaryFilePath("ramd_groups.dat");
    TextWriter::writeFileFromString(groupsFilePath, c_groupsFileContent);

    KeyValueTreeBuilder mdpValueBuilder;
    mdpValueBuilder.rootObject().addValue(std::string(RAMDModuleInfo::sc_name) + "-active",
                                          std::string("true"));
    mdpValueBuilder.rootObject().addValue(std::string(RAMDModuleInfo::sc_name) + "-groups-file",
                                          groupsFilePath.string());
    RAMDOptions ramdOptions;
    test::fillOptionsFromMdpValues(mdpValueBuilder.build(), &ramdOptions);
    // No mdp directory set, matching the pre-existing (cwd-relative/absolute) behavior
    EXPECT_NO_THROW(ramdOptions.setInputGroupIndices(ramdIndexGroupsAndNames()));
    EXPECT_EQ(1, ramdOptions.parameters().ngroups_);
}

} // namespace
} // namespace gmx
