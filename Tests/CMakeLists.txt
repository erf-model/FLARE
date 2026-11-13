# FLARE's tests: the regression tests of Tests/CTestList.cmake (FLARE_ENABLE_TESTS) and the
# GoogleTest unit tests of Tests/Unit (FLARE_ENABLE_UNIT_TESTS). They are registered on
# FLARE's own executables, built from ERF's library; ERF's own tests stay in ERF.

set(ERF_TEST_NRANKS "2" CACHE STRING "Number of MPI ranks to use for each regression test")
set(ERF_TEST_FCOMPARE_RTOL "2.0e-10" CACHE STRING "fcompare relative tolerance")
set(ERF_TEST_FCOMPARE_ATOL "2.0e-10" CACHE STRING "fcompare absolute tolerance")
set(ERF_TEST_GOLD_FILES_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/ERFGoldFiles" CACHE PATH "Path to FLARE's gold files")

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
  set(ERF_TEST_PYTHON "${Python3_EXECUTABLE}")
else()
  set(ERF_TEST_PYTHON "")
endif()

# ERF's general helpers (resolve_test_exe, setup_test, add_test_r, add_test_abort,
# add_test_restart_abort); each runs ${ERF_TEST_EXE}, here FLARE's executable
include(${ERF_SOURCE_DIR}/Tests/ERFTestHelpers.cmake)
set(ERF_TEST_EXE flare)

if(FLARE_ENABLE_TESTS)
  # ERF is added EXCLUDE_FROM_ALL, so the AMReX plotfile tools the tests run are built with
  # FLARE's default target
  add_custom_target(flare_test_tools ALL)
  foreach(_tool IN ITEMS fcompare fextract fextrema fvolumesum)
    if(TARGET ${_tool})
      add_dependencies(flare_test_tools ${_tool})
    endif()
  endforeach()

  # flare_examples: ERF's example model interfaces, each with one example model -- the fire
  # engine (FireEngineExample), the dust model (DustModelExample) and ERF's structure model
  # (Submodules/ERF/Tests/StructureModelExample) -- and FLARE's problem file. It is separate
  # from flare, so that flare keeps one fire engine and one dust model.
  set(FLARE_EXAMPLE_FIRE ${ERF_ENABLE_FIRE})
  set(FLARE_EXAMPLE_DUST ${ERF_ENABLE_DUST})
  set(FLARE_EXAMPLE_STRUCTURES ${ERF_ENABLE_MOVING_BODIES})
  if(FLARE_EXAMPLE_FIRE OR FLARE_EXAMPLE_DUST OR FLARE_EXAMPLE_STRUCTURES)
    add_executable(flare_examples)
    target_sources(flare_examples PRIVATE ${PROJECT_SOURCE_DIR}/Exec/ERF_Prob.cpp)
    target_include_directories(flare_examples PRIVATE ${PROJECT_SOURCE_DIR}/Exec)
    if(FLARE_EXAMPLE_FIRE)
      set(ERF_FIRE_ENGINE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/FireEngineExample")
      include("${ERF_FIRE_ENGINE_DIR}/ERFFireEngine.cmake")
      erf_fire_engine_add(flare_examples)
    endif()
    if(FLARE_EXAMPLE_DUST)
      set(ERF_DUST_MODEL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/DustModelExample")
      include("${ERF_DUST_MODEL_DIR}/ERFDustModel.cmake")
      erf_dust_model_add(flare_examples)
    endif()
    if(FLARE_EXAMPLE_STRUCTURES)
      set(ERF_STRUCTURE_MODEL_DIR "${ERF_SOURCE_DIR}/Tests/StructureModelExample")
      include("${ERF_STRUCTURE_MODEL_DIR}/ERFStructureModel.cmake")
      erf_structure_model_add(flare_examples)
    endif()
    erf_setup_executable(flare_examples)
  endif()

  include(${CMAKE_CURRENT_SOURCE_DIR}/CTestList.cmake)
endif()

if(FLARE_ENABLE_UNIT_TESTS)
  add_subdirectory(Unit)
endif()
