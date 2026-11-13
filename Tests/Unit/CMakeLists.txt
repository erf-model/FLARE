# FLARE's GoogleTest unit tests: the fire engine (Fire/), the dust layer (Dust/) and the
# structures (Structures/Tests/Unit), in FLARE's own executable, flare_unit_tests, with
# FLARE's models compiled in (flare_add_models) and ERF's test main and shared helpers
# (Submodules/ERF/Tests/Unit). ERF's own unit tests stay in ERF.

set(BUILD_GMOCK OFF CACHE BOOL "Build Google Mock targets for the vendored GoogleTest dependency")
set(INSTALL_GTEST OFF CACHE BOOL "Install targets from the vendored GoogleTest dependency")
if(MSVC)
  set(gtest_force_shared_crt ON CACHE BOOL "Use the shared CRT when building the vendored GoogleTest dependency")
endif()
if(NOT TARGET GTest::gtest)
  add_subdirectory(${ERF_SOURCE_DIR}/Submodules/googletest
                   ${PROJECT_BINARY_DIR}/Submodules/googletest
                   EXCLUDE_FROM_ALL)
endif()
include(GoogleTest)

set(FLARE_GTEST_REPORT_ROOT "${PROJECT_BINARY_DIR}/test-results/gtest")
set(_erf_unit_dir ${ERF_SOURCE_DIR}/Tests/Unit)

set(flare_unit_exe flare_unit_tests)
add_executable(${flare_unit_exe})
target_sources(${flare_unit_exe} PRIVATE
  ${PROJECT_SOURCE_DIR}/Exec/ERF_Prob.cpp
  ${_erf_unit_dir}/ERF_GTestMain.cpp
  ${_erf_unit_dir}/ERF_GTestRuntimeStubs.cpp
)
target_include_directories(${flare_unit_exe} PRIVATE
  ${PROJECT_SOURCE_DIR}/Exec
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${_erf_unit_dir}
  ${_erf_unit_dir}/MovingBodies)

# The fire headers are only on the include path, and the fire sources only in
# the library, when the fire module is built.
if(ERF_ENABLE_FIRE)
  target_sources(${flare_unit_exe} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestLcpReader.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestFarsiteWindOutput.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestFarsiteSpreadAccumulation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestFireAcceleration.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestFireBoundaryGuard.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestFireAnchorLevel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestLevelSetEllipse.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestBehaveTransfer.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestBehaveScottBurgan.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestFbp.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestFuelMoisture.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestFuelMoistureStick.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestFuelMoistureEMC.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestFuelMoistureLive.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestFirePrecip.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestLevelSetAdvection.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestPrescribedFire.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestDirectionalShape.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestSplitHamiltonianAdvection.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestWindLimit.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestReactionVelocityFormula.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestFuelMap.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestStructureIgnition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestFireSuppression.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestCustomFuel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestFireDiagnosticsPerFuel.cpp
    # Motivation: several perimeter files in one deck must all be read and merged.
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestPolygonIgnition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestFireModelReferences.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestFireDefaults.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestFireModelSummary.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestBehaveNetLoad.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestCheneyGould.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestFarsiteShape.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestFireHeatBudget.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestArrivalStamp.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestFbpSlopeVector.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestFireBoxParity.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestFireAtmInjection.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestFireCoupler.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestFireRegistry.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestRothermelCellMoisture.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fire/ERF_GTestEmberDrift.cpp
  )
endif()

# The dust column kernels (settling, deposition, emission): headers on the
# include path only when the dust module is built.
if(ERF_ENABLE_DUST)
  target_sources(${flare_unit_exe} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/Dust/ERF_GTestDustColumn.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Dust/ERF_GTestDustThreshold.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Dust/ERF_GTestDustWind.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Dust/ERF_GTestDustBudget.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Dust/ERF_GTestDustDiagnostics.cpp
  )
endif()

if(ERF_ENABLE_STRUCTURES)
  # FLARE's structures (Structures/): the MoorDyn wrapper and the tower models
  include(${PROJECT_SOURCE_DIR}/Structures/Tests/Unit/ERFStructuresUnitTests.cmake)
  erf_structures_unit_tests(${flare_unit_exe})
endif()

# FLARE's models, then ERF's library, flags and (with NO_ERF_MAIN) no ERF main.cpp
flare_add_models(${flare_unit_exe})
erf_setup_executable(${flare_unit_exe} NO_ERF_MAIN)
target_link_libraries(${flare_unit_exe} PRIVATE GTest::gtest)

file(MAKE_DIRECTORY "${FLARE_GTEST_REPORT_ROOT}/${flare_unit_exe}")
gtest_discover_tests(${flare_unit_exe}
  DISCOVERY_MODE POST_BUILD
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DISCOVERY_TIMEOUT 60
  XML_OUTPUT_DIR "${FLARE_GTEST_REPORT_ROOT}/${flare_unit_exe}"
  PROPERTIES
    LABELS unit
    TIMEOUT 600
    SKIP_RETURN_CODE 77
)

add_test(NAME flare_unit_tests_shuffle_repeat
  COMMAND ${flare_unit_exe} --gtest_shuffle --gtest_repeat=3
          --gtest_output=xml:${FLARE_GTEST_REPORT_ROOT}/${flare_unit_exe}/shuffle-repeat.xml)
set_tests_properties(flare_unit_tests_shuffle_repeat PROPERTIES LABELS unit TIMEOUT 1800)

# The fire and dust kernels under the floating-point traps the fire decks run
# with. clang at -O2/-O3 evaluates the unselected arm of `c ? a/b : x`, so a
# guarded 0/0 or log(0) raises FE_INVALID/FE_DIVBYZERO and kills the run; the
# plain unit run never arms the traps and cannot see it (before the October
# 2026 audit BehaveTransfer.WindowReachesState died here with SIGILL in the
# BEHAVE live-weight 0/0 of a fuel without live load).
if(ERF_ENABLE_FIRE)
  set(_fire_gtest_suites "BehaveScottBurgan.*:BehaveTransfer.*:CustomFuel.*:DirectionalShape.*:FarsiteFrontCell.*:FarsiteSpreadTest.*:Fbp.*:FireAcceleration.*:FireAnchorLevel.*:FireBoundaryGuard.*:FireDiagnosticsPerFuel.*:FirePrecip.*:FireSuppression.*:FuelMap.*:FuelMoisture.*:FuelMoistureEMC.*:FuelMoistureLive.*:FuelMoistureStick.*:LcpReader.*:LevelSetAdvection.*:LevelSetEllipse.*:PrescribedFire.*:SplitHamiltonianAdvection.*:StructureIgnition.*:WindLimit.*:FireModelReferences.*:FireDefaults.*:BehaveNetLoad.*:CheneyGould98.*:MacArthurCap.*:FarsiteShape.*:FireHeatBudget.*:ByramIntensity.*:ArrivalStamp.*:FbpSlopeVector.*:PerimeterCount.*:FuelBlendingParity.*:CrownFront.*:ExnerTendency.*:SmokeInjection.*:RothermelCellMoisture.*:EmberDrift.*:FireBreak.*:ReactionVelocityFormula.*:WrfBmstCompat.*:FireModelSummary.*:FireCoupler.*:FireRegistry.*:HeatInjection.*")
  if(ERF_ENABLE_DUST)
    string(APPEND _fire_gtest_suites ":DustColumn.*:DustThreshold.*:DustWind.*:DustBudget.*:DustDiagnostics.*")
  endif()
  add_test(NAME FireDustGTests_FpeTraps
    COMMAND ${flare_unit_exe}
      --gtest_filter=${_fire_gtest_suites}
      amrex.fpe_trap_invalid=1 amrex.fpe_trap_zero=1 amrex.fpe_trap_overflow=1)
  set_tests_properties(FireDustGTests_FpeTraps PROPERTIES LABELS "unit;fire" TIMEOUT 600)
endif()

# The structures' one collective, the SubDyn reader (the I/O rank reads, every rank gets the
# broadcast), on 2 ranks, with ERF's parallel test main
if(ERF_ENABLE_STRUCTURES AND ERF_ENABLE_MPI)
  set(flare_parallel_exe flare_parallel_tests)
  add_executable(${flare_parallel_exe})
  target_sources(${flare_parallel_exe} PRIVATE
    ${PROJECT_SOURCE_DIR}/Exec/ERF_Prob.cpp
    ${_erf_unit_dir}/ERF_GTestParallelMain.cpp
    ${_erf_unit_dir}/ERF_GTestRuntimeStubs.cpp)
  target_include_directories(${flare_parallel_exe} PRIVATE ${PROJECT_SOURCE_DIR}/Exec ${_erf_unit_dir})
  erf_structures_parallel_tests(${flare_parallel_exe})
  flare_add_models(${flare_parallel_exe})
  erf_setup_executable(${flare_parallel_exe} NO_ERF_MAIN)
  target_link_libraries(${flare_parallel_exe} PRIVATE GTest::gtest)
  add_test(NAME flare_parallel_tests_np2
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${flare_parallel_exe}>)
  set_tests_properties(flare_parallel_tests_np2 PROPERTIES LABELS "unit;parallel" PROCESSORS 2 TIMEOUT 600)
endif()
