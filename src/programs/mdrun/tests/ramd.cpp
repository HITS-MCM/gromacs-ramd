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

TEST_F(RAMDTest, RAMDCanRun)
{
    runner_.useTopGroAndNdxFromDatabase("spc2");
    const std::string mdpContents = R"(
        dt                       = 0.004
        nsteps                   = 2
        tcoupl                   = Berendsen
        tc-grps                  = System
        tau-t                    = 0.5
        ref-t                    = 300
        constraints              = all-bonds
        cutoff-scheme            = Verlet
    )";
//        ramd                     = yes
//        ramd-seed                = 9876
//        ramd-force               = 600.0
//        ramd-eval-freq           = 50
//        ramd-r-min-dist          = 0.025
//        ramd-force-out-freq      = 100
//        ramd-max-dist            = 4.0
//        ramd-receptor            = FirstWaterMolecule
//        ramd-ligand              = SecondWaterMolecule
//    )";
    runner_.useStringAsMdpFile(mdpContents);

    EXPECT_EQ(0, runner_.callGrompp());

    ::gmx::test::CommandLine caller;

    // Do an mdrun with RAMD enabled
    ASSERT_EQ(0, runner_.callMdrun(caller));
}

} // namespace test
} // namespace gmx
