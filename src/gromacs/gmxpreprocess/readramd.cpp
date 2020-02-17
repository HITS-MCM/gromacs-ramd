#include "gmxpre.h"

#include <cassert>
#include <cstdlib>
#include <cstring>

#include "gromacs/domdec/localatomsetmanager.h"
#include "gromacs/fileio/readinp.h"
#include "gromacs/fileio/warninp.h"
#include "gromacs/gmxpreprocess/readir.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/ramd_params.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

char **read_ramdparams(std::vector<t_inpfile> *inp,
                       gmx::RAMDParams *ramdparams,
                       warninp_t wi)
{
    ramdparams->seed = get_eint(inp, "ramd-seed", 1234, wi);
    ramdparams->force = get_ereal(inp, "ramd-force", 600, wi);
	ramdparams->eval_freq = get_eint(inp, "ramd-eval-freq", 50, wi);
	ramdparams->r_min_dist = get_ereal(inp, "ramd-r-min-dist", 0.0025, wi);
	ramdparams->force_out_freq = get_eint(inp, "ramd-force-out-freq", 100, wi);
	ramdparams->max_dist = get_ereal(inp, "ramd-max-dist", 4.0, wi);

    char buf[STRLEN];
    char **grpbuf;
    snew(grpbuf, 2);
    snew(grpbuf[0], STRLEN);
    snew(grpbuf[1], STRLEN);

    sprintf(buf, "ramd-receptor");
    setStringEntry(inp, buf, grpbuf[0], "");
    sprintf(buf, "ramd-ligand");
    setStringEntry(inp, buf, grpbuf[1], "");

	ramdparams->old_angle_dist = (get_eeenum(inp, "ramd-old-angle-dist", yesno_names, wi) != 0);

    return grpbuf;
}

void make_ramd_groups(gmx::RAMDParams *ramdparams,
                      char **pgnames,
                      const t_blocka *grps, char **gnames)
{
    int ig = search_string(pgnames[0], grps->nr, gnames);
    ramdparams->protein.nat = grps->index[ig + 1] - grps->index[ig];

    snew(ramdparams->protein.ind, ramdparams->protein.nat);
	for (int i = 0; i < ramdparams->protein.nat; i++)
	{
		ramdparams->protein.ind[i] = grps->a[grps->index[ig] + i];
    }

    ig = search_string(pgnames[1], grps->nr, gnames);
    ramdparams->ligand.nat = grps->index[ig + 1] - grps->index[ig];

    snew(ramdparams->ligand.ind, ramdparams->ligand.nat);
	for (int i = 0; i < ramdparams->ligand.nat; i++)
	{
		ramdparams->ligand.ind[i] = grps->a[grps->index[ig] + i];
    }
}
