# The fragility fit (Structures/Tools/fragility_fit/fragility_fit.py): from the demand table an
# ensemble gathers (Structures/Tools/fragility_ensemble) it fits a lognormal fragility function
# per component and per group of components by maximum likelihood (Baker 2015, multiple stripes
# analysis, doi:10.1193/021113EQS025M). Tests/CTestList.cmake includes this file when
# ERF_ENABLE_STRUCTURES is ON. The tests are Python only (standard library) and run no ERF
# executable, hence the "unit" label.
#   FragilityFit_Published  the standard deviations and coefficients of variation the paper prints
#                           for its simulated studies (p. 589 and Table 3), and its Eq. 1
#   FragilityFit_ByHand     fits that follow from Phi(0) and Phi(1), a brute-force search of the
#                           likelihood, the outcomes without a fit, a lognormal capacity
#   FragilityFit_Table      two stored demand tables through the command: components, groups,
#                           directions, the points per stripe
#   FragilityFit_Refused    every plan key out of range and a table that cannot be used stop the
#                           tool, named
#   FragilityFit_DamageStates  DS1-DS4 counted by hand from a damage table, DS1 equal to the group
#                           of every component, the cascade threshold, a member that skips a
#                           state, and every damage_states key out of range
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
  foreach(_ff_class Published ByHand Table Refused DamageStates)
    add_test(NAME FragilityFit_${_ff_class}
             COMMAND "${Python3_EXECUTABLE}" -I
                     "${PROJECT_SOURCE_DIR}/Structures/Tools/fragility_fit/tests/test_fragility_fit.py"
                     ${_ff_class})
    set_tests_properties(FragilityFit_${_ff_class}
        PROPERTIES
        TIMEOUT 300
        PROCESSORS 1
        LABELS "unit;conductors;fragility")
  endforeach()
else()
  message(STATUS "FragilityFit_* need Python 3 and are not registered")
endif()
