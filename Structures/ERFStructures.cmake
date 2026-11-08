# FLARE's structures (Structures/): the MoorDyn-C wrapper and the lattice tower models.
# The top-level CMakeLists.txt includes this file when ERF_ENABLE_STRUCTURES is ON.
# It turns on the actuator core and the structure coupling (ERF_ENABLE_MOVING_BODIES),
# adds this directory to ERF_STRUCTURE_MODEL_DIRS, so ERF compiles the sources listed in
# ERFStructureModel.cmake, and finds MoorDyn-C or builds the bundled stub.

list(APPEND CMAKE_MESSAGE_CONTEXT "Structures")
message(STATUS "Configuring FLARE structures (Structures/)")

# The structures sit on the actuator core and plug into ERF's structure coupling: turn
# ERF_ENABLE_MOVING_BODIES on for this configure (a normal variable shadowing the cache entry,
# so a later configure with ERF_ENABLE_STRUCTURES=OFF gets the cached choice back)
if(NOT ERF_ENABLE_MOVING_BODIES)
  message(STATUS "ERF_ENABLE_STRUCTURES=ON turns ERF_ENABLE_MOVING_BODIES on for this configure")
endif()
set(ERF_ENABLE_MOVING_BODIES ON)
# Structures/ is compiled in through this option only: a user entry in
# ERF_STRUCTURE_MODEL_DIRS that resolves to it (relative path, trailing slash) is refused
# rather than added twice
foreach(_erf_dir IN LISTS ERF_STRUCTURE_MODEL_DIRS)
  get_filename_component(_erf_dir_abs "${_erf_dir}" ABSOLUTE)
  if(_erf_dir_abs STREQUAL "${PROJECT_SOURCE_DIR}/Structures")
    message(FATAL_ERROR "ERF_STRUCTURE_MODEL_DIRS lists ${_erf_dir}, FLARE's Structures/ directory; "
                        "ERF_ENABLE_STRUCTURES=ON compiles it in, so leave it out of ERF_STRUCTURE_MODEL_DIRS")
  endif()
endforeach()
list(APPEND ERF_STRUCTURE_MODEL_DIRS "${PROJECT_SOURCE_DIR}/Structures")

set(MOORDYN_DIR "" CACHE PATH "MoorDyn-C install prefix (Structures/setup_moordyn.sh); read unless ERF_MOORDYN_USE_STUB=ON")

if(ERF_MOORDYN_USE_STUB)
  # A stand-in for libmoordyn implementing the subset of the MoorDyn-C v2 C API that the
  # wrapper calls, so the structures build and their tests run where MoorDyn-C is not
  # installed; it has no line dynamics (see MoorDyn/Stub/ERF_MoorDynStub.cpp).
  message(STATUS "Using the bundled MoorDyn stub library (ERF_MOORDYN_USE_STUB=ON)")
  add_library(erf_moordyn_stub STATIC ${PROJECT_SOURCE_DIR}/Structures/MoorDyn/Stub/ERF_MoorDynStub.cpp)
  target_include_directories(erf_moordyn_stub PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/Structures/MoorDyn/Stub>
    $<INSTALL_INTERFACE:include>)
  # position independent, as erf_srclib links it PUBLIC and may itself be shared (BUILD_SHARED_LIBS)
  set_target_properties(erf_moordyn_stub PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON
                                                    POSITION_INDEPENDENT_CODE ON)
  if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)   # installed with ERF, as CMakeLists.txt does
    install(TARGETS erf_moordyn_stub EXPORT ERFTargets ARCHIVE DESTINATION lib)
    # the installed ERF_MoorDynSystem.H includes <moordyn/MoorDyn2.h>: install the stub's copy with it
    install(FILES ${PROJECT_SOURCE_DIR}/Structures/MoorDyn/Stub/moordyn/MoorDyn2.h
            DESTINATION include/moordyn)
  endif()
  set(ERF_MOORDYN_VERSION "stub")
else()
  # MoorDyn-C installs lib/cmake/moordyn/MoorDynConfig.cmake (target MoorDyn::moordyn) and
  # its headers under include/moordyn; MOORDYN_DIR is the install prefix
  # (Structures/setup_moordyn.sh).
  find_package(MoorDyn 2.3 REQUIRED CONFIG HINTS ${MOORDYN_DIR})
  message(STATUS "Found MoorDyn-C ${MoorDyn_VERSION}: ${MoorDyn_DIR}")
  # The wind reaches the lines through the v2 external wave-kinematics calls; make sure the
  # installed header declares them (they replaced the v1 names in 2.3).
  get_target_property(_erf_moordyn_incs MoorDyn::moordyn INTERFACE_INCLUDE_DIRECTORIES)
  # looked up afresh each configure, so a changed MOORDYN_DIR is not checked against an old header
  unset(ERF_MOORDYN_C_HEADER CACHE)
  find_file(ERF_MOORDYN_C_HEADER MoorDyn2.h
    PATHS ${_erf_moordyn_incs} ${MOORDYN_DIR}/include
    PATH_SUFFIXES moordyn
    NO_DEFAULT_PATH)
  if(NOT ERF_MOORDYN_C_HEADER)
    message(FATAL_ERROR "moordyn/MoorDyn2.h not found under ${_erf_moordyn_incs}")
  endif()
  file(STRINGS "${ERF_MOORDYN_C_HEADER}" _erf_moordyn_extkin REGEX "MoorDyn_ExternalWaveKinSet")
  if(NOT _erf_moordyn_extkin)
    message(FATAL_ERROR "FLARE's structures need MoorDyn-C 2.3 or newer: MoorDyn_ExternalWaveKinSet "
                        "is not declared in ${ERF_MOORDYN_C_HEADER}")
  endif()
  message(STATUS "MoorDyn C API header: ${ERF_MOORDYN_C_HEADER}")
  set(ERF_MOORDYN_VERSION "${MoorDyn_VERSION}")
endif()

message(STATUS "FLARE structures configuration complete")
list(POP_BACK CMAKE_MESSAGE_CONTEXT)
