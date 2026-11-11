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
  # Motivation: the conductor inputs must be read with their defaults and refuse every value
  # outside its documented range with a message naming the key; the MoorDyn input written from
  # them must carry the span lowered by the surface offset, the air density, external kinematics
  # and fixed attachments, and MoorDyn must accept it; a span must hang at the catenary sag in
  # still air, report positions in ERF's frame, and blow out towards atan(q/w) in a prescribed
  # crosswind, stepping in lockstep with ERF's clock; the manager must place the attachments
  # at their height above the terrain under each end, step on the anchor level only, and hand
  # MoorDyn ERF's velocity sampled where the line is now. A section over towers must be read,
  # written with its insulator strings and free points, and hang from the strings; the closest
  # approach of two lines must be the exact minimum distance between their conductors, flagged
  # against the flashover distance. The line ends on a transformer's footprint must be dead-ended
  # on it, and its load must be the sum of the lines' pulls and their moment about its base,
  # flagged against its allowable values; the closest approach of a conductor to its box must be
  # the exact distance. The drag's support must hold every spread point with its own reach.
  # Each component with a capacity must be checked once per anchor-level step (never on another
  # level), its demand pinned to the still-air catenary tension or the footings' share of the
  # weight, its first exceedance logged once and kept across a restart, and failure_stop must stop
  # the run after the log is written (ERF_GTestCapacities.cpp holds the monitor's own tests).
  # With intensity_measure = wind the averaged wind beside each component's peak must be the time-weighted mean
  # over the averaging window, kept at the peak's step (not the last), and continue across a restart to the bit
  # (ERF_GTestIntensity.cpp).
  # With ground_contact the terrain under each line must reach MoorDyn as its seafloor, node for
  # node, and hold a slack span on the ground; the contact's damping must be critical for the
  # conductor and a MoorDyn step too long for it refused; the spans' ground records must survive a
  # restart (ERF_GTestGroundContact.cpp holds the log's own tests).
  # With failure_break a line's component over its capacity must let go of the line at the end of that step, at
  # the end the rule names, once; the released span must come to hang from its other end with the weight of its
  # hanging length there, its support take nothing, and a restart make the release again and continue the run
  # (ERF_GTestLineBreaks.cpp holds the rule's, the log's and the line's own tests).
  # With tower_failure = release a tower over a footing's or member's capacity must let go of every line on it, once,
  # at the end of that step, carry no line afterwards and stay failed across a restart; on MoorDyn-C the next tower,
  # which then carries more of the conductor, must fail after it (ERF_GTestTowerFailures.cpp holds the rule's and
  # the log's own tests).
  # After a release the section's state at rest must be the one the balance of forces gives: the residual tension in
  # the spans that stay, the lean of the strings and the towers' pulls, against an independent solver
  # (Structures/Tools/broken_wire; ERF_GTestBrokenWire.cpp, MoorDyn-C only).
  target_sources(${target} PRIVATE
    ${_erf_structures_tests_dir}/ERF_GTestConductorInputsReset.cpp
    ${_erf_structures_tests_dir}/ERF_GTestConductorGeometry.cpp
    ${_erf_structures_tests_dir}/ERF_GTestConductorInputs.cpp
    ${_erf_structures_tests_dir}/ERF_GTestASCE74.cpp
    ${_erf_structures_tests_dir}/ERF_GTestGusts.cpp
    ${_erf_structures_tests_dir}/ERF_GTestMoorDynInputWriter.cpp
    ${_erf_structures_tests_dir}/ERF_GTestConductorLine.cpp
    ${_erf_structures_tests_dir}/ERF_GTestConductors.cpp
    ${_erf_structures_tests_dir}/ERF_GTestCapacities.cpp
    ${_erf_structures_tests_dir}/ERF_GTestIntensity.cpp
    ${_erf_structures_tests_dir}/ERF_GTestGroundContact.cpp
    ${_erf_structures_tests_dir}/ERF_GTestLineBreaks.cpp
    ${_erf_structures_tests_dir}/ERF_GTestTowerFailures.cpp
    ${_erf_structures_tests_dir}/ERF_GTestBrokenWire.cpp
    ${_erf_structures_tests_dir}/ERF_GTestTransformer.cpp
  )
  # Motivation: a conductor's temperature under a fire (conductors.thermal_model = ieee738) must use the
  # IEEE 738-2023 terms as an independent Python calculation writes them out
  # (Tests/test_files/Conductors_FireHeating/check_fire_heating.py), its point-source irradiance must give
  # the closed-form line integral over a long conductor, every new key must be range-checked, and the span's
  # length must follow its mean temperature in MoorDyn and come back from a checkpoint.
  target_sources(${target} PRIVATE
    ${_erf_structures_tests_dir}/ERF_GTestLineThermal.cpp
    ${_erf_structures_tests_dir}/ERF_GTestTowerThermal.cpp
  )
  # Motivation: against the real MoorDyn-C (the tests skip on the stub, which has no line
  # dynamics), a span's mean swing in a steady crosswind must be the quasi-static blowout angle
  # atan(q/w), its free swing after a gust must have the cable's first out-of-plane period
  # 2 c / sqrt(H/m), and the suspension strings of a section must swing to the angle of the wind
  # span over the weight span, each within 1 %; the pulls on a level span's dead ends must carry
  # its weight and the elastic catenary's horizontal tension.
  target_sources(${target} PRIVATE
    ${_erf_structures_tests_dir}/ERF_GTestConductorVerification.cpp
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
