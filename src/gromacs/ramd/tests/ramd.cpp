/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017 by the GROMACS development team.
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include "gmxpre.h"

#include "gromacs/ramd/ramd.h"

#include <gtest/gtest.h>

#include "gromacs/commandline/filenm.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pulling/pull_internal.h"

namespace gmx
{
namespace test
{
namespace
{

struct RAMDTest : public ::testing::Test
{
protected:

    RAMDParams params;
    std::unique_ptr<RAMD> ramd;
    CommrecHandle cr;

    void SetUp() override
    {
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

        cr = init_commrec(MPI_COMM_WORLD, nullptr);

        t_filenm fnm[] = {
            { efXVG, "-ramd", "ramd", ffOPTWR }
        };

        int64_t step = 0;
        pull_t pull;
        t_pull_coord coord_params;
        coord_params.eType = epullEXTERNAL;
        char buf[] = "RAMD";
        coord_params.externalPotentialProvider = buf;
        pull_coord_work_t pull_coord(coord_params);
        pull.coord.push_back(pull_coord);
        pull.coord.push_back(pull_coord);
        pull.coord.push_back(pull_coord);
        pull.numUnregisteredExternalPotentials = 3;
        pull.ePBC = 0;

        ramd = std::make_unique<RAMD>(params, &pull, &step, StartingBehavior::NewSimulation, cr.get(), 1, fnm, nullptr);
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
    t_mdatoms md;
    double time = 0.0;
    matrix box = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    ForceProviderInput forceProviderInput(x, md, time, box, *cr);

    PaddedVector<RVec> f = {{0, 0, 0}};
    ForceWithVirial forceWithVirial(f, true);
    gmx_enerdata_t enerd(1, 0);
    ForceProviderOutput forceProviderOutput(&forceWithVirial, &enerd);

    ramd->calculateForces(forceProviderInput, &forceProviderOutput);
}

} // namespace
} // namespace test
} // namespace gmx
