/**
 * \internal \file
 *
 * \brief
 * Implements test of RAMD
 *
 * \author Bernd Doser <bernd.doser@h-its.org>
 *
 */

#include "gmxpre.h"
#include "programs/mdrun/tests/moduletest.h"

namespace gmx {
namespace test {

class RAMDTestFixture : public MdrunTestFixture
{
    protected:
        RAMDTestFixture();
        ~RAMDTestFixture() override;
};


RAMDTestFixture::RAMDTestFixture()
{}

RAMDTestFixture::~RAMDTestFixture()
{}


//! Test fixture for mdrun with RAMD settings
typedef gmx::test::RAMDTestFixture RAMDTest;

TEST_F(RAMDTest, RAMD_1whhi)
{
    runner_.useTopGroAndNdxFromDatabase("1wdhi");
    const std::string mdpContents = R"(
        dt                       = 0.004
        nsteps                   = 2
        tcoupl                   = Berendsen
        tc-grps                  = System
        tau-t                    = 0.5
        ref-t                    = 300
        constraints              = all-bonds
        cutoff-scheme            = Verlet
        ramd                     = yes
        ramd-seed                = 9876
        ramd-ngroups             = 1
        ramd-group1-receptor     = Protein
        ramd-group1-ligand       = INH
        ramd-group1-force        = 600.0
        ramd-group1-r-min-dist   = 0.025
        ramd-group1-max-dist     = 4.0
        ramd-eval-freq           = 50
        ramd-force-out-freq      = 0
    )";
    runner_.useStringAsMdpFile(mdpContents);

    EXPECT_EQ(0, runner_.callGrompp());

    // Do an mdrun with RAMD enabled
    ASSERT_EQ(0, runner_.callMdrun());
}

TEST_F(RAMDTest, RAMD_membrane)
{
    runner_.useTopGroAndNdxFromDatabase("membrane");
    const std::string mdpContents = R"(
        dt                       = 0.004
        nsteps                   = 2
        tcoupl                   = Berendsen
        tc-grps                  = System
        tau-t                    = 0.5
        ref-t                    = 300
        constraints              = all-bonds
        cutoff-scheme            = Verlet
        ramd                     = yes
        ramd-seed                = 9876
        ramd-ngroups             = 1
        ramd-group1-receptor     = Protein
        ramd-group1-ligand       = IXO
        ramd-group1-force        = 600.0
        ramd-group1-r-min-dist   = 0.025
        ramd-group1-max-dist     = 4.0
        ramd-eval-freq           = 50
        ramd-force-out-freq      = 0
    )";
    runner_.useStringAsMdpFile(mdpContents);

    EXPECT_EQ(0, runner_.callGrompp());

    // Do an mdrun with RAMD enabled
    ASSERT_EQ(0, runner_.callMdrun());
}

// heat shock protein 90
TEST_F(RAMDTest, RAMD_hsp90)
{
    runner_.useTopGroAndNdxFromDatabase("hsp90");
    const std::string mdpContents = R"(
        integrator               = md
        dt                       = 0.002
        nsteps                   = 200
        nstlog                   = 2500
        nstenergy                = 2500
        nstxout-compressed       = 500
        continuation             = yes
        constraints              = h-bonds
        constraint-algorithm     = lincs
        cutoff-scheme            = Verlet
        coulombtype              = PME
        rcoulomb                 = 1.2
        vdwtype                  = Cut-off
        DispCorr                 = EnerPres
        ramd                     = yes
        ramd-seed                = 989
        ramd-ngroups             = 1
        ramd-group1-receptor     = Protein
        ramd-group1-ligand       = INH
        ramd-group1-force        = 585.0
        ramd-group1-r-min-dist   = 0.0025
        ramd-group1-max-dist     = 4.0
        ramd-eval-freq           = 60
        ramd-force-out-freq      = 0
        comm-mode                = Linear
        nstcomm                  = 100
        comm_grps                = System
        tinit                    = 0.000
        nstxout                  = 25000
        nstvout                  = 50000
        compressed-x-precision   = 1000
        compressed-x-grps        = SYSTEM
        pbc                      = xyz
        rlist                    = 1.10
        fourierspacing           = 0.12
        pme_order                = 4
        ewald_geometry           = 3d
        ewald-rtol               = 1e-5
        ewald-rtol-lj            = 1e-5
        vdw-modifier             = Potential-shift
        rvdw-switch              = 0.0
        rvdw                     = 1.2
        tcoupl                   = Nose-Hoover
        tc_grps                  = Protein_INH Water_and_ions
        tau_t                    = 1.0  1.0
        ref_t                    = 300  300
        Pcoupl                   = Parrinello-Rahman
        pcoupltype               = semiisotropic
        tau_p                    = 5.0
        compressibility          = 4.5e-5  4.5e-5
        ref_p                    = 1.0     1.0
        gen_vel                  = no
    )";
    runner_.useStringAsMdpFile(mdpContents);

    EXPECT_EQ(0, runner_.callGrompp());

    // Do an mdrun with RAMD enabled
    ASSERT_EQ(0, runner_.callMdrun());
}

} // namespace test
} // namespace gmx
