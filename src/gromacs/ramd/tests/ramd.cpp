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
#include "gromacs/gmxpreprocess/readir.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pulling/pull_internal.h"
#include "gromacs/ramd/ramd.h"
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

struct RAMDTest : public ::testing::Test
{
protected:

    std::unique_ptr<RAMD> ramd;
    std::unique_ptr<t_commrec> cr = std::make_unique<t_commrec>();

    void SetUp() override
    {
        // t_inputrec ir;
        // ir.bPull = true;
        // ir.pull = std::make_unique<pull_params_t>();

        // t_pull_coord coord_params;
        // coord_params.eType = PullingAlgorithm::External;
        // coord_params.externalPotentialProvider = "RAMD";
        // ir.pull->coord.push_back(coord_params);
        // ir.pull->coord.push_back(coord_params);
        // ir.pull->coord.push_back(coord_params);

        // t_pull_group pull_group;
        // ir.pull->ngroup = 1;
        // ir.pull->group.push_back(pull_group);

        // gmx_mtop_t mtop;
        // addNWaterMolecules(&mtop, 2);
        // mtop.finalize();

        // t_state state;

        // CommandLine cmdline;
        // cmdline.addOption("grompp");

        gmx_mtop_t mtop;
        t_inputrec ir;
        t_state state;
        read_tpx_state(TestFileManager::getInputFilePath("data/1WDHI/topol.tpr").u8string(), &ir, &state, &mtop);

        WarningHandler wi{true, 0};

        pull_t* pull = set_pull_init(
            &ir, mtop, state.x, state.box, state.lambda[FreeEnergyPerturbationCouplingType::Mass], &wi);

        int64_t step = 0;
        RAMDParams params;
        params.seed = 9876;
        params.ngroup = 1;
        snew(params.group, 1);
        params.group->force = 600.0;
        params.group->max_dist = 4.0;
        params.group->r_min_dist = 0.025;
        params.eval_freq = 50;
        params.force_out_freq = 50;
        params.old_angle_dist = false;
        params.eval_freq = 1;

        t_filenm fnm[] = {
            { efXVG, "-ramd", "ramd", ffOPTWR }
        };

        ramd = std::make_unique<RAMD>(params, pull, &step, StartingBehavior::NewSimulation,
            cr.get(), 1, fnm, nullptr);
    }

    void TearDown() override
    {
        ramd.reset();
    }
};

TEST_F(RAMDTest, construction)
{
    EXPECT_TRUE(ramd);
}

TEST_F(RAMDTest, calculateForces)
{
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
