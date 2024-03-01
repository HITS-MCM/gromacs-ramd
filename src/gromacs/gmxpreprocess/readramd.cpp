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

#include "gromacs/fileio/readinp.h"
#include "gromacs/fileio/warninp.h"
#include "gromacs/gmxpreprocess/readir.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/ramd_params.h"
#include "gromacs/mdtypes/pull_params.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/smalloc.h"

void read_ramdparams(std::vector<t_inpfile>* inp, gmx::RAMDParams* ramdparams, WarningHandler* wi)
{
    ramdparams->seed = get_eint(inp, "ramd-seed", 1234, wi);
    ramdparams->ngroup = get_eint(inp, "ramd-ngroups", 1, wi);

    for (int i = 0; i < ramdparams->ngroup; i++)
    {
        gmx::RAMDGroup ramdgrp;
        auto ramd_prefix = std::string("ramd-group") + std::to_string(i + 1);
        auto pull1_prefix = std::string("pull-group") + std::to_string(i * 2 + 1);
        auto pull2_prefix = std::string("pull-group") + std::to_string(i * 2 + 2);

        inp->emplace_back(0, 1, false, false, false, pull1_prefix + "-name",
            get_estr(inp, ramd_prefix + "-receptor", "Protein"));
        inp->emplace_back(0, 1, false, false, false, pull1_prefix + "-pbcatom",
            get_estr(inp, ramd_prefix + "-receptor-pbc-atom", "0"));

        inp->emplace_back(0, 1, false, false, false, pull2_prefix + "-name",
            get_estr(inp, ramd_prefix + "-ligand", "INH"));
        inp->emplace_back(0, 1, false, false, false, pull2_prefix + "-pbcatom",
            get_estr(inp, ramd_prefix + "-ligand-pbc-atom", "0"));

        ramdgrp.force = get_ereal(inp, ramd_prefix + "-force", 600, wi);
        ramdgrp.r_min_dist = get_ereal(inp, ramd_prefix + "-r-min-dist", 0.0025, wi);
        ramdgrp.max_dist = get_ereal(inp, ramd_prefix + "-max-dist", 4.0, wi);
        ramdgrp.bind_res_receptor = get_estr(inp, ramd_prefix + "-receptor-res", "bind_res_receptor");
        ramdgrp.bind_res_ligand = get_estr(inp, ramd_prefix + "-ligand-res", "bind_res_ligand");

        ramdparams->group.emplace_back(ramdgrp);
    }

    inp->emplace_back(0, 1, false, false, false, "pull-pbc-ref-prev-step-com",
        get_estr(inp, "ramd-pbc-ref-prev-step-com", "yes"));

    ramdparams->eval_freq = get_eint(inp, "ramd-eval-freq", 50, wi);
    ramdparams->force_out_freq = get_eint(inp, "ramd-force-out-freq", 100, wi);
    ramdparams->old_angle_dist = getEnum<Boolean>(inp, "ramd-old-angle-dist", wi) != Boolean::No;

    inp->emplace_back(0, 1, false, false, false, "pull-ngroups",
        std::to_string(ramdparams->ngroup * 2));
    inp->emplace_back(0, 1, false, false, false, "pull-nstxout",
        std::to_string(ramdparams->force_out_freq));
    inp->emplace_back(0, 1, false, false, false, "pull-nstfout",
        std::to_string(ramdparams->force_out_freq));
    inp->emplace_back(0, 1, false, false, false, "pull-ncoords",
        std::to_string(ramdparams->ngroup * 3));

    std::vector<std::string> v{"1 0 0", "0 1 0", "0 0 1"};
    for (int i = 0; i < ramdparams->ngroup; ++i)
    {
        for (int d = 0; d < 3; ++d)
        {
            auto prefix = std::string("pull-coord") + std::to_string(i * 3 + d + 1);
            inp->emplace_back(0, 1, false, false, false,
                prefix + "-groups", std::to_string(i * 2 + 1) + " " + std::to_string(i * 2 + 2));
            inp->emplace_back(0, 1, false, false, false,
                prefix + "-type", "external-potential");
            inp->emplace_back(0, 1, false, false, false,
                prefix + "-potential-provider", "RAMD");
            inp->emplace_back(0, 1, false, false, false,
                prefix + "-geometry", "direction");
            inp->emplace_back(0, 1, false, false, false,
                prefix + "-vec", v[d]);
        }
    }

    ramdparams->connected_ligands = getEnum<Boolean>(inp, "ramd-connected-ligands", wi) != Boolean::No;
    ramdparams->residence_distance = get_ereal(inp, "residence-distance", 0.55, wi);
}


void add_residence_time_groups(t_inputrec* ir, std::vector<IndexGroup> indexGroups)
{
    if (find_group(ir->ramdParams->group[0].bind_res_receptor.c_str(), indexGroups) == -1) return;

    const int gid = getGroupIndex(ir->ramdParams->group[0].bind_res_receptor, indexGroups);
    auto group = indexGroups[gid].particleIndices;

    t_pull_group new_group;
    new_group.ind.push_back(0);
    new_group.ind.push_back(1);
    new_group.ind.push_back(2);
    new_group.pbcatom = 1;
    ir->pull->group.push_back(new_group);
    ir->pull->ngroup++;

    t_pull_group new_group2;
    new_group2.ind.push_back(3);
    new_group2.ind.push_back(4);
    new_group2.ind.push_back(5);
    new_group2.pbcatom = 4;
    ir->pull->group.push_back(new_group2);
    ir->pull->ngroup++;

    t_pull_coord new_coord;
    new_coord.ngroup = 2;
    new_coord.group[0] = ir->ramdParams->ngroup * 2 + 1;
    new_coord.group[1] = ir->ramdParams->ngroup * 2 + 2;
    new_coord.coordIndex = ir->pull->ncoord;
    new_coord.dim = {1, 1, 1};
    new_coord.eGeom = PullGroupGeometry::Distance;
    ir->pull->coord.push_back(new_coord);
    ir->pull->ncoord++;
}
