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

#include "gromacs/fileio/readinp.h"
#include "gromacs/fileio/warninp.h"
#include "gromacs/gmxpreprocess/readir.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/ramd_params.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/smalloc.h"

char** read_ramdparams(std::vector<t_inpfile>* inp, gmx::RAMDParams* ramdparams, warninp_t wi)
{
    ramdparams->seed = get_eint(inp, "ramd-seed", 1234, wi);
    ramdparams->ngroup = get_eint(inp, "ramd-ngroups", 1, wi);

    char buf[STRLEN];
    char** grpbuf;
    snew(grpbuf, ramdparams->ngroup);
    for (int i = 1; i < ramdparams->ngroup; i++)
    {
        auto ramdgrp = &ramdparams->group[i];
        snew(grpbuf[i], STRLEN);
        sprintf(buf, "ramd-group%d-receptor", i);
        setStringEntry(inp, buf, grpbuf[i], "");
        sprintf(buf, "ramd-group%d-ligand", i);
        setStringEntry(inp, buf, grpbuf[i], "");
        sprintf(buf, "ramd-group%d-force", i);
        ramdgrp->force = get_ereal(inp, buf, 600, wi);
        sprintf(buf, "ramd-group%d-r-min-dist", i);
        ramdgrp->r_min_dist = get_ereal(inp, buf, 0.0025, wi);
        sprintf(buf, "ramd-group%d-max-dist", i);
        ramdgrp->max_dist  = get_ereal(inp, buf, 4.0, wi);
    }

    ramdparams->eval_freq = get_eint(inp, "ramd-eval-freq", 50, wi);
    ramdparams->force_out_freq = get_eint(inp, "ramd-force-out-freq", 100, wi);
    ramdparams->old_angle_dist = (get_eeenum(inp, "ramd-old-angle-dist", yesno_names, wi) != 0);

    snew(grpbuf, 2);
    snew(grpbuf[0], STRLEN);
    snew(grpbuf[1], STRLEN);

    sprintf(buf, "ramd-receptor");
    setStringEntry(inp, buf, grpbuf[0], "");
    sprintf(buf, "ramd-ligand");
    setStringEntry(inp, buf, grpbuf[1], "");

    return grpbuf;
}

void make_ramd_groups(gmx::RAMDParams* ramdparams, char** pgnames, const t_blocka* grps, char** gnames)
{
    int ig                  = search_string(pgnames[0], grps->nr, gnames);
    ramdparams->protein.nat = grps->index[ig + 1] - grps->index[ig];

    snew(ramdparams->protein.ind, ramdparams->protein.nat);
    for (int i = 0; i < ramdparams->protein.nat; i++)
    {
        ramdparams->protein.ind[i] = grps->a[grps->index[ig] + i];
    }

    ig                     = search_string(pgnames[1], grps->nr, gnames);
    ramdparams->ligand.nat = grps->index[ig + 1] - grps->index[ig];

    snew(ramdparams->ligand.ind, ramdparams->ligand.nat);
    for (int i = 0; i < ramdparams->ligand.nat; i++)
    {
        ramdparams->ligand.ind[i] = grps->a[grps->index[ig] + i];
    }
}
