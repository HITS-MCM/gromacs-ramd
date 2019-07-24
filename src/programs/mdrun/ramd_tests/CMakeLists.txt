set(testname "RAMDIntegrationTests")
set(exename "ramd-integration-tests")

gmx_add_gtest_executable(
    ${exename}
    # files with code for tests
    RAMDTest.cpp
    # pseudo-library for code for mdrun
    $<TARGET_OBJECTS:mdrun_objlib>
)
    
gmx_register_gtest_test(
    ${testname}
    ${exename}
    INTEGRATION_TEST
)
