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
        ramd-force               = 600.0
        ramd-eval-freq           = 50
        ramd-r-min-dist          = 0.025
        ramd-force-out-freq      = 100
        ramd-max-dist            = 4.0
        ramd-receptor            = Protein
        ramd-ligand              = INH
    )";
    runner_.useStringAsMdpFile(mdpContents);

    EXPECT_EQ(0, runner_.callGrompp());

    ::gmx::test::CommandLine caller;

    // Do an mdrun with RAMD enabled
    ASSERT_EQ(0, runner_.callMdrun(caller));
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
        ramd-force               = 600.0
        ramd-eval-freq           = 50
        ramd-r-min-dist          = 0.025
        ramd-force-out-freq      = 100
        ramd-max-dist            = 4.0
        ramd-receptor            = Protein
        ramd-ligand              = IXO
    )";
    runner_.useStringAsMdpFile(mdpContents);

    EXPECT_EQ(0, runner_.callGrompp());

    ::gmx::test::CommandLine caller;

    // Do an mdrun with RAMD enabled
    ASSERT_EQ(0, runner_.callMdrun(caller));
}

} // namespace test
} // namespace gmx
