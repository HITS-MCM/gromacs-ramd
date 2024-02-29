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
#include "gmxpre.h"

#include <gtest/gtest.h>

#include "gromacs/commandline/filenm.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/warninp.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/gmxpreprocess/readir.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pulling/pull_internal.h"
#include "gromacs/ramd/ramd.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testfilemanager.h"
#include "testutils/topologyhelpers.h"

namespace gmx
{
namespace test
{
namespace
{

TEST(RAMDTest, CalculateForces1WDHI)
{
    gmx_mtop_t mtop;
    t_inputrec ir;
    t_state state;
    read_tpx_state(TestFileManager::getInputFilePath("data/1WDHI/topol.tpr").u8string(), &ir, &state, &mtop);

    WarningHandler wi{true, 0};

    pull_t* pull = set_pull_init(
        &ir, mtop, state.x, state.box, state.lambda[FreeEnergyPerturbationCouplingType::Mass], &wi);

    t_filenm fnm[] = {
        { efXVG, "-ramd", "ramd", ffOPTWR }
    };

    auto cr = std::make_unique<t_commrec>();
    auto ramd = std::make_unique<RAMD>(*ir.ramdParams, pull, StartingBehavior::NewSimulation,
        cr.get(), 1, fnm, nullptr, mtop);

    ASSERT_NEAR(0.0, pull->coord[0].scalarForce, 1e-6);

    PaddedVector<RVec> x = {{0, 0, 0}};
    std::vector<real> chargeA{1};
    matrix box = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    ForceProviderInput forceProviderInput(x, ssize(chargeA), chargeA, {}, 0.0, 0, box, *cr);

    PaddedVector<RVec> f = {{0, 0, 0}};
    ForceWithVirial forceWithVirial(f, true);
    gmx_enerdata_t enerd(1, 0);
    ForceProviderOutput forceProviderOutput(&forceWithVirial, &enerd);

    ramd->calculateForces(forceProviderInput, &forceProviderOutput);

    ASSERT_NEAR(528.87046337127686, pull->coord[0].scalarForce, 1e-6);
}

TEST(RAMDTest, binding_residues)
{
    CommandLine cmdline;
    cmdline.addOption("grompp");
    cmdline.addOption("-f", TestFileManager::getInputFilePath("data/4water.mdp").u8string());
    cmdline.addOption("-c", TestFileManager::getInputFilePath("data/4water.gro").u8string());
    cmdline.addOption("-p", TestFileManager::getInputFilePath("data/4water.top").u8string());
    cmdline.addOption("-n", TestFileManager::getInputFilePath("data/4water.ndx").u8string());
    gmx::test::TestFileManager fileManager;
    std::string outTprFilename = fileManager.getTemporaryFilePath("4water.tpr").u8string();
    cmdline.addOption("-o", outTprFilename);

    ASSERT_EQ(0, gmx_grompp(cmdline.argc(), cmdline.argv()));

    gmx_mtop_t mtop;
    t_inputrec ir;
    t_state state;
    read_tpx_state(outTprFilename.c_str(), &ir, &state, &mtop);

    WarningHandler wi{true, 0};

    pull_t* pull = set_pull_init(
        &ir, mtop, state.x, state.box, state.lambda[FreeEnergyPerturbationCouplingType::Mass], &wi);

    t_filenm fnm[] = {
        { efXVG, "-ramd", "ramd", ffOPTWR }
    };

    auto cr = std::make_unique<t_commrec>();
    auto ramd = std::make_unique<RAMD>(*ir.ramdParams, pull, StartingBehavior::NewSimulation,
        cr.get(), 1, fnm, nullptr, mtop);

    ASSERT_EQ(7, pull->params.group.size());
    ASSERT_EQ(0, pull->params.group[0].ind.size());

    ASSERT_EQ((std::vector<int>{0, 1, 2}), pull->params.group[1].ind);
    ASSERT_EQ((std::vector<int>{3, 4, 5}), pull->params.group[2].ind);
    ASSERT_EQ((std::vector<int>{0, 1, 2}), pull->params.group[3].ind);
    ASSERT_EQ((std::vector<int>{6, 7, 8}), pull->params.group[4].ind);

    ASSERT_EQ("1SOL", ramd->getParams().group[0].bind_res_receptor);
    ASSERT_EQ("2SOL", ramd->getParams().group[0].bind_res_ligand);

    t_atoms atoms = gmx_mtop_global_atoms(mtop);
    ASSERT_EQ(atoms.atom[0].resind, 0);
    ASSERT_EQ(atoms.nres, 4);
    ASSERT_EQ(std::string(atoms.resinfo[0].name[0]), "SOL");
    ASSERT_EQ(std::string(atoms.resinfo[1].name[0]), "SOL");
    ASSERT_EQ(std::string(atoms.atomtype[0][0]), "OW_spc");

    PaddedVector<RVec> x = {{0, 0, 0}};
    std::vector<real> chargeA{1};
    matrix box = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    ForceProviderInput forceProviderInput(x, ssize(chargeA), chargeA, {}, 0.0, 0, box, *cr);

    PaddedVector<RVec> f = {{0, 0, 0}};
    ForceWithVirial forceWithVirial(f, true);
    gmx_enerdata_t enerd(1, 0);
    ForceProviderOutput forceProviderOutput(&forceWithVirial, &enerd);

    ramd->calculateForces(forceProviderInput, &forceProviderOutput);
}

} // namespace
} // namespace test
} // namespace gmx
