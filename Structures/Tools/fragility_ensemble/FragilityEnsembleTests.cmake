# The fragility ensemble driver (Structures/Tools/fragility_ensemble/fragility_ensemble.py): from a
# base deck and a plan of wind-speed stripes, wind directions and gust seeds it writes one run
# directory per member, runs the members behind a load gate and gathers their demand and intensity
# tables (conductors.intensity_measure = wind). Tests/CTestList.cmake includes this file when
# ERF_ENABLE_STRUCTURES is ON. The tests are Python only (standard library) and run no ERF
# executable, hence the "unit" label; their expected decks, wind files, manifest and collected
# tables are written out by hand in tests/test_fragility_ensemble.py.
#   FragilityEnsemble_DryRun   the member set, each member's deck and wind files, and the manifest
#                              of a dry run on the stored base decks (inflow through a face,
#                              periodic, prescribed wind)
#   FragilityEnsemble_Refused  every plan and run option out of range stops the driver, named
#   FragilityEnsemble_Run      a stand-in executable: finished members are skipped, failed ones run
#                              again, the load gate waits, pauses and resumes
#   FragilityEnsemble_Collect  stored member tables, with a member without a table, a failed one,
#                              one not run and one with another header
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND AND NOT WIN32)
  foreach(_fe_class DryRun Refused Run Collect)
    add_test(NAME FragilityEnsemble_${_fe_class}
             COMMAND "${Python3_EXECUTABLE}" -I
                     ${PROJECT_SOURCE_DIR}/Structures/Tools/fragility_ensemble/tests/test_fragility_ensemble.py
                     ${_fe_class})
    set_tests_properties(FragilityEnsemble_${_fe_class}
        PROPERTIES
        TIMEOUT 120
        PROCESSORS 1
        LABELS "unit;conductors;fragility")
  endforeach()
else()
  message(STATUS "FragilityEnsemble_* need Python 3 and a POSIX system and are not registered")
endif()
