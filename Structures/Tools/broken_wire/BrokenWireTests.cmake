# The state at rest of a conductor section after a release
# (Structures/Tools/broken_wire/broken_wire.py): every span and string a heavy elastic cable, the
# forces at every string's lower end in balance. It is the reference, independent of ERF and
# MoorDyn, that the unit tests BrokenWire.* compare MoorDyn-C at rest with.
# Tests/CTestList.cmake includes this file when ERF_ENABLE_STRUCTURES is ON. The tests are Python
# only (standard library) and run no ERF executable, hence the "unit" label.
#   BrokenWire_Cable      a cable's reach against a numerical integration of its shape, the lumped
#                         chain against the cable, the level span against the elastic catenary of
#                         the documentation
#   BrokenWire_Balance    solved sections against checks that do not use the tool's balance at
#                         the strings' feet: the supports carry the weight, the potential energy
#                         is stationary, a string points along its load
#   BrokenWire_Limits     states known without the tool
#   BrokenWire_Robust     a section whose strings lean far from plumb; a conductor that would lie on
#                         the ground is refused
#   BrokenWire_Reference  the numbers ERF_GTestBrokenWire.cpp is pinned to, and the command
#   BrokenWire_Refused    every plan key out of range stops the tool, named
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
  foreach(_bw_class Cable Balance Limits Robust Reference Refused)
    add_test(NAME BrokenWire_${_bw_class}
             COMMAND "${Python3_EXECUTABLE}" -I
                     "${PROJECT_SOURCE_DIR}/Structures/Tools/broken_wire/tests/test_broken_wire.py"
                     ${_bw_class})
    set_tests_properties(BrokenWire_${_bw_class}
        PROPERTIES
        TIMEOUT 300
        PROCESSORS 1
        LABELS "unit;conductors;failure")
  endforeach()
else()
  message(STATUS "BrokenWire_* need Python 3 and are not registered")
endif()
