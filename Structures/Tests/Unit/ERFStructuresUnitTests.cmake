# The unit tests of FLARE's structures (Structures/), added to ERF's gtest executables by
# Tests/Unit/CMakeLists.txt when ERF_ENABLE_STRUCTURES is ON. The tests include the shared
# helpers in Tests/Unit (ERF_GTestTempDir.H) and Tests/Unit/MovingBodies
# (ERF_GTestThrowOnAbort.H).
set(_erf_structures_tests_dir ${CMAKE_CURRENT_LIST_DIR})
set(_erf_structures_frame_files ${CMAKE_CURRENT_LIST_DIR}/../test_files/FrameSubDynTower)

function(erf_structures_unit_tests target)
  target_include_directories(${target} PRIVATE
    ${_erf_structures_tests_dir}
    ${PROJECT_SOURCE_DIR}/Tests/Unit
    ${PROJECT_SOURCE_DIR}/Tests/Unit/MovingBodies)
  # Motivation: a lattice tower's member drag must be the normal-flow drag of every node, equal
  # the hand values in a uniform wind and the fine integral in a log-law wind; the one-mode tower
  # must respond to a step load with its frequency and damping and its foundation must take the
  # drag nodes' loads; every tower-type input outside its range must be refused by name.
  target_sources(${target} PRIVATE
    ${_erf_structures_tests_dir}/ERF_GTestTowers.cpp
    ${_erf_structures_tests_dir}/ERF_GTestTowerDynamics.cpp
  )
  # Motivation: a lattice tower's frame model must reproduce the closed-form static response of beams
  # (cantilever bending with and without shear deformation, axial and torsional stiffness, a
  # spring base, a fixed portal frame, a triangulated truss), balance every load with its support
  # reactions, refuse a mechanism naming the free degree of freedom, read a SubDyn input file and
  # refuse what it cannot model, and give the stiffness at a lattice tower's peak that SubDyn
  # gives for the same input file (Structures/Tests/test_files/FrameSubDynTower). Its mass and modes
  # must match beam theory and SubDyn's natural frequencies, and its Newmark response the exact
  # discrete response of each mode, without losing energy when undamped. A tower standing on its
  # frame must carry its drag and line loads onto the frame without losing force or moment, settle
  # where the frame's static solution puts it, and give its footings the loads and its weight.
  target_sources(${target} PRIVATE
    ${_erf_structures_tests_dir}/ERF_GTestFrame.cpp
    ${_erf_structures_tests_dir}/ERF_GTestFrameSubDyn.cpp
    ${_erf_structures_tests_dir}/ERF_GTestFrameDynamics.cpp
    ${_erf_structures_tests_dir}/ERF_GTestFrameTower.cpp
    ${_erf_structures_tests_dir}/ERF_GTestMemberChecks.cpp
  )
  # Motivation: the MoorDyn system wrapper must report the configured version and a missing input
  # file, initialise a fixed-fixed span with no coupled degree of freedom at the catenary sag and
  # tension, expose the external kinematics points starting with the line nodes, blow the span out
  # under a crosswind towards the quasi-static angle, continue identically from a saved state, and
  # run under the floating-point traps. The bundled stub or the installed MoorDyn-C answers,
  # whichever the build links.
  target_sources(${target} PRIVATE
    ${_erf_structures_tests_dir}/ERF_GTestMoorDynSystem.cpp
  )
  target_compile_definitions(${target} PRIVATE ERF_FRAME_TEST_FILES="${_erf_structures_frame_files}")
endfunction()

function(erf_structures_parallel_tests target)
  target_include_directories(${target} PRIVATE
    ${_erf_structures_tests_dir}
    ${PROJECT_SOURCE_DIR}/Tests/Unit
    ${PROJECT_SOURCE_DIR}/Tests/Unit/MovingBodies)
  # Motivation: the SubDyn file reader is the tower code's one collective (the I/O rank reads,
  # every rank receives the broadcast): every rank must get the same frame and the same refusals.
  target_sources(${target} PRIVATE
    ${_erf_structures_tests_dir}/ERF_GTestFrameSubDyn.cpp
  )
  target_compile_definitions(${target} PRIVATE ERF_FRAME_TEST_FILES="${_erf_structures_frame_files}")
endfunction()
