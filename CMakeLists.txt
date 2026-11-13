############################ BASE ######################################

# FLARE uses ERF the way ERF uses AMReX: ERF (Submodules/ERF) is a library, ERF::erf, and
# FLARE builds its own executable, flare, from its own sources -- the fire engine
# (FireEngine/), the dust layer (Dust/), the structures (Structures/) and its problem file
# (Exec/ERF_Prob.cpp) -- with ERF's main.cpp, compile flags and GPU settings
# (erf_setup_executable). FLARE hands ERF none of its model-directory hooks
# (ERF_FIRE_ENGINE_DIRS, ERF_DUST_MODEL_DIRS, ERF_STRUCTURE_MODEL_DIRS): its models register
# themselves (ERF_REGISTER_FIRE_ENGINE, _DUST_MODEL, _STRUCTURE_MODEL) from files compiled
# into flare.
#
#   git submodule update --init --recursive
#   cmake -S . -B build -DERF_ENABLE_MPI=ON [-DERF_ENABLE_DUST=ON] [-DERF_ENABLE_STRUCTURES=ON]
#   cmake --build build
#
# Every ERF option (ERF_ENABLE_MPI, ERF_ENABLE_CUDA, ERF_ENABLE_NETCDF, ...) is given on this
# configure line as for ERF itself.

cmake_minimum_required(VERSION 3.25 FATAL_ERROR)
project(FLARE CXX C)

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/Submodules/ERF/CMake/BuildERFExe.cmake")
  message(FATAL_ERROR "Submodules/ERF is empty: run git submodule update --init --recursive")
endif()

########################## OPTIONS #####################################

# FLARE's defaults for ERF's options: the fire coupling is on, since FLARE is a fire model.
# These are cache entries made before ERF's own option() calls, so a value on the configure
# line still wins.
option(ERF_ENABLE_FIRE "Enable the fire coupling and FLARE's fire engine (FireEngine/)" ON)
option(ERF_ENABLE_DUST "Enable the dust coupling and FLARE's dust layer (Dust/)" OFF)

# FLARE's own options
option(ERF_ENABLE_STRUCTURES "Build FLARE's structures in Structures/ (the MoorDyn-C wrapper, the tower models and the conductor model; implies ERF_ENABLE_MOVING_BODIES)" OFF)
option(ERF_MOORDYN_USE_STUB "Build FLARE's structures against the bundled MoorDyn stub instead of MoorDyn-C (for testing where MoorDyn-C is not installed)" OFF)
if(ERF_MOORDYN_USE_STUB AND NOT ERF_ENABLE_STRUCTURES)
  message(FATAL_ERROR "ERF_MOORDYN_USE_STUB needs ERF_ENABLE_STRUCTURES=ON: the stub stands in for MoorDyn-C in FLARE's structures")
endif()
option(FLARE_ENABLE_TESTS "Enable FLARE's regression tests (Tests/CTestList.cmake)" OFF)
option(FLARE_ENABLE_UNIT_TESTS "Enable FLARE's GoogleTest unit tests (Tests/Unit)" ${FLARE_ENABLE_TESTS})

# ERF's own tests stay in ERF: FLARE registers only its own
set(ERF_ENABLE_TESTS OFF)
set(ERF_ENABLE_UNIT_TESTS OFF)
if(FLARE_ENABLE_TESTS)
  set(ERF_ENABLE_FCOMPARE ON)   # amrex_fcompare and the other plotfile tools FLARE's tests use
endif()

if(ERF_ENABLE_FIRE)
  # nvlink refused the fire kernels and the line and plane samplers with "uses too much data
  # for compiler-generated constants" (October 2026). ERF's set_erf_compile_flags appends to
  # this list, for ERF's library and for flare alike.
  list(APPEND ERF_CUDA_FLAGS "--ptxas-options=--disable-optimizer-constants")
endif()

########################## STRUCTURES ##################################

# MoorDyn-C (or the stub) and ERF_ENABLE_MOVING_BODIES, before ERF is configured
if(ERF_ENABLE_STRUCTURES)
  include(${PROJECT_SOURCE_DIR}/Structures/ERFStructures.cmake)
endif()

############################ ERF #######################################

# CUDA is enabled here as well as in ERF: enable_language in ERF's directory does not reach
# this one, where flare and the unit-test executable are defined, and CMake then stops with
# "No known features for CUDA compiler".
if(ERF_ENABLE_CUDA)
  enable_language(CUDA)
endif()

enable_testing()
add_subdirectory(Submodules/ERF EXCLUDE_FROM_ALL)

########################## EXECUTABLE ##################################

set(flare_exe_name flare)
add_executable(${flare_exe_name})
target_sources(${flare_exe_name} PRIVATE ${PROJECT_SOURCE_DIR}/Exec/ERF_Prob.cpp)
target_include_directories(${flare_exe_name} PRIVATE ${PROJECT_SOURCE_DIR}/Exec)

# FLARE's models, compiled into a target. Each directory's file defines the function that
# adds its sources to a target; *_DIR names the directory during the call. They are added
# before erf_setup_executable, which compiles every .cpp of the target as CUDA in a CUDA
# build. Tests/Unit calls this for FLARE's unit-test executable too.
function(flare_add_models target)
  if(ERF_ENABLE_FIRE)
    set(ERF_FIRE_ENGINE_DIR "${PROJECT_SOURCE_DIR}/FireEngine")
    include("${ERF_FIRE_ENGINE_DIR}/ERFFireEngine.cmake")
    erf_fire_engine_add(${target})
  endif()
  if(ERF_ENABLE_DUST)
    set(ERF_DUST_MODEL_DIR "${PROJECT_SOURCE_DIR}/Dust")
    include("${ERF_DUST_MODEL_DIR}/ERFDustModel.cmake")
    erf_dust_model_add(${target})
  endif()
  if(ERF_ENABLE_STRUCTURES)
    set(ERF_STRUCTURE_MODEL_DIR "${PROJECT_SOURCE_DIR}/Structures")
    include("${ERF_STRUCTURE_MODEL_DIR}/ERFStructureModel.cmake")
    erf_structure_model_add(${target})
  endif()
endfunction()

flare_add_models(${flare_exe_name})
erf_setup_executable(${flare_exe_name})

############################ TESTS #####################################

if(FLARE_ENABLE_TESTS OR FLARE_ENABLE_UNIT_TESTS)
  add_subdirectory(Tests)
endif()
