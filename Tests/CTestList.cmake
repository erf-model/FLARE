# FLARE's regression tests. Tests/CMakeLists.txt sets the test variables and includes
# ERF's general helpers (Submodules/ERF/Tests/ERFTestHelpers.cmake: resolve_test_exe,
# setup_test, add_test_r, add_test_abort, add_test_restart_abort) before this file, with
# ERF_TEST_EXE = flare, so every helper below runs FLARE's executable. ERF's own tests stay
# in ERF; the files FLARE's tests take from ERF (a start-up deck, the plotfile reader, the
# shared RANS checks, ERF's parity drivers) are read from ${ERF_SOURCE_DIR}.

#=============================================================================
# Helpers copied from ERF's Tests/CTestList.cmake, which FLARE's conductor tests use
#=============================================================================

# Run the deck in test_files/<TEST_FILES_DIR> twice, once with a diagnostic off and once
# with it on, and require PLTFILE to be identical bit for bit.  This is the harness for
# "turning this output on does not change the answer": OFF_OPTIONS and ON_OPTIONS are the
# two option sets, COMMON_OPTIONS go to both, and REQUIRE_ON_FILE names a file the on leg
# must write and the off leg must not, so a misspelled option cannot pass as agreement.
function(add_test_option_parity TEST_NAME TEST_FILES_DIR PLTFILE)
    set(oneValueArgs "COMMON_OPTIONS" "OFF_OPTIONS" "ON_OPTIONS" "REQUIRE_ON_FILE" "RUN_TIMEOUT")
    cmake_parse_arguments(ADD_TEST_OP "" "${oneValueArgs}" "" ${ARGN})
    setup_test()
    resolve_test_exe("" "${ERF_TEST_EXE}" TEST_EXE)

    set(_run_timeout 600)
    set(_ctest_timeout 600)
    if(DEFINED ADD_TEST_OP_RUN_TIMEOUT)
        set(_run_timeout "${ADD_TEST_OP_RUN_TIMEOUT}")
        math(EXPR _ctest_timeout "2 * ${_run_timeout} + 600")
    endif()

    add_test(NAME ${TEST_NAME} COMMAND ${CMAKE_COMMAND}
        "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
        "-DNRANKS=${NP}"
        "-DTEST_EXE=${TEST_EXE}"
        "-DCONFIG=$<CONFIG>"
        "-DINPUT=${CURRENT_TEST_BINARY_DIR}/${TEST_FILES_DIR}.i"
        "-DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}"
        "-DFCOMPARE=${FCOMPARE_EXE}"
        "-DPLTFILE=${PLTFILE}"
        "-DRUN_TIMEOUT=${_run_timeout}"
        "-DCOMMON_OPTIONS=${ADD_TEST_OP_COMMON_OPTIONS}"
        "-DOFF_OPTIONS=${ADD_TEST_OP_OFF_OPTIONS}"
        "-DON_OPTIONS=${ADD_TEST_OP_ON_OPTIONS}"
        "-DREQUIRE_ON_FILE=${ADD_TEST_OP_REQUIRE_ON_FILE}"
        -P ${ERF_SOURCE_DIR}/Tests/RunOptionParity.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT ${_ctest_timeout}
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;option-parity"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/option_off/simulation.log;${CURRENT_TEST_BINARY_DIR}/option_on/simulation.log;${CURRENT_TEST_BINARY_DIR}/parity.log")
endfunction(add_test_option_parity)

# Restart parity: run one deck straight, then to a checkpoint and on from it, and
# require the plotfile at the end to be identical (no gold file). Every run has a
# time limit; the default stays at 600, but an explicit RUN_TIMEOUT is forwarded
# unchanged to each leg and used to size the outer CTest watchdog.
# COMMON_OPTIONS reaches all three legs; RESTART_OPTIONS reaches the restart leg alone,
# which is where anything that changes the decomposition has to go.  ALLOW_DIFF_GRIDS
# lets fcompare compare plotfiles written on different BoxArrays, which a restart leg
# that re-makes the level-0 grids needs and no other restart test should want.
# CHK_NRANKS/RESTART_NRANKS give the checkpoint and restart legs separate rank counts, for
# the case of continuing a run on more ranks than it was written with; the straight leg
# follows the restart's count so the comparison isolates the restart, not the decomposition.
# REQUIRE_LEVEL0_REMAKE makes a regrid test fail if the restart stopped regridding, which
# would otherwise leave it silently comparing an ordinary restart and passing.
function(add_test_restart_parity TEST_NAME TEST_FILES_DIR STEP_CHK STEP_END)
    set(oneValueArgs "COMMON_OPTIONS" "RESTART_OPTIONS" "CHK_OPTIONS" "CHK_LEG_END" "CHK_NRANKS" "RESTART_NRANKS" "FCOMPARE_RTOL" "FCOMPARE_ATOL" "RUN_TIMEOUT" "DATALOG" "DATALOG_SIGDIGITS" "PLT2DFILE" "BNDRY_PLANES_DIR" "BNDRY_PLANES_FIRST_STEP" "BNDRY_PLANES_STEPS" "BNDRY_PLANES_READ_INPUT" "CHECKER" "CHECK_MODE" "CHECKER_LOG")
    cmake_parse_arguments(ADD_TEST_RP "ALLOW_DIFF_GRIDS;REQUIRE_LEVEL0_REMAKE;BNDRY_PLANES_STALE;BNDRY_PLANES_CUT_NEWLINE;BNDRY_PLANES_EARLY_STOP" "${oneValueArgs}" "" ${ARGN})
    setup_test()
    resolve_test_exe("" "${ERF_TEST_EXE}" TEST_EXE)
    set(_checker_log "")
    set(_post_checkpoint_log "")
    if(NOT "${ADD_TEST_RP_CHECK_MODE}" STREQUAL "")
        if("${ADD_TEST_RP_CHECKER}" STREQUAL "")
            message(FATAL_ERROR "${TEST_NAME}: CHECKER is required when CHECK_MODE is set")
        endif()
        if("${ADD_TEST_RP_CHECKER_LOG}" STREQUAL "")
            set(_checker_log "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.checkpoint-activity.log")
        else()
            set(_checker_log "${ADD_TEST_RP_CHECKER_LOG}")
        endif()
        if(ADD_TEST_RP_CHECK_MODE STREQUAL "wdm6_active")
            set(_post_checkpoint_log "${CURRENT_TEST_BINARY_DIR}/post_checkpoint_parity.log")
        endif()
    endif()

    set(_fcompare_rtol "${ERF_TEST_FCOMPARE_RTOL}")
    set(_fcompare_atol "${ERF_TEST_FCOMPARE_ATOL}")
    if(NOT "${ADD_TEST_RP_FCOMPARE_RTOL}" STREQUAL "")
        set(_fcompare_rtol "${ADD_TEST_RP_FCOMPARE_RTOL}")
    endif()
    if(NOT "${ADD_TEST_RP_FCOMPARE_ATOL}" STREQUAL "")
        set(_fcompare_atol "${ADD_TEST_RP_FCOMPARE_ATOL}")
    endif()
    set(_run_timeout 600)
    set(_ctest_timeout 600)
    if(DEFINED ADD_TEST_RP_RUN_TIMEOUT)
        set(_run_timeout "${ADD_TEST_RP_RUN_TIMEOUT}")
        math(EXPR _ctest_timeout "3 * ${_run_timeout} + 600")
    endif()

    add_test(NAME ${TEST_NAME} COMMAND ${CMAKE_COMMAND}
        "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
        "-DNRANKS=${NP}"
        "-DTEST_EXE=${TEST_EXE}"
        "-DCONFIG=$<CONFIG>"
        "-DINPUT=${CURRENT_TEST_BINARY_DIR}/${TEST_FILES_DIR}.i"
        "-DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}"
        "-DFCOMPARE=${FCOMPARE_EXE}"
        "-DSTEP_CHK=${STEP_CHK}"
        "-DSTEP_END=${STEP_END}"
        "-DRTOL=${_fcompare_rtol}"
        "-DATOL=${_fcompare_atol}"
        "-DRUN_TIMEOUT=${_run_timeout}"
        "-DCOMMON_OPTIONS=${ADD_TEST_RP_COMMON_OPTIONS}"
        "-DRESTART_OPTIONS=${ADD_TEST_RP_RESTART_OPTIONS}"
        "-DALLOW_DIFF_GRIDS=${ADD_TEST_RP_ALLOW_DIFF_GRIDS}"
        "-DREQUIRE_LEVEL0_REMAKE=${ADD_TEST_RP_REQUIRE_LEVEL0_REMAKE}"
        "-DCHK_NRANKS=${ADD_TEST_RP_CHK_NRANKS}"
        "-DRESTART_NRANKS=${ADD_TEST_RP_RESTART_NRANKS}"
        "-DDATALOG=${ADD_TEST_RP_DATALOG}"
        "-DDATALOG_SIGDIGITS=${ADD_TEST_RP_DATALOG_SIGDIGITS}"
        "-DPLT2DFILE=${ADD_TEST_RP_PLT2DFILE}"
        "-DCHK_OPTIONS=${ADD_TEST_RP_CHK_OPTIONS}"
        "-DCHK_LEG_END=${ADD_TEST_RP_CHK_LEG_END}"
        "-DBNDRY_PLANES_DIR=${ADD_TEST_RP_BNDRY_PLANES_DIR}"
        "-DBNDRY_PLANES_FIRST_STEP=${ADD_TEST_RP_BNDRY_PLANES_FIRST_STEP}"
        "-DBNDRY_PLANES_STALE=${ADD_TEST_RP_BNDRY_PLANES_STALE}"
        "-DBNDRY_PLANES_STEPS=${ADD_TEST_RP_BNDRY_PLANES_STEPS}"
        "-DBNDRY_PLANES_CUT_NEWLINE=${ADD_TEST_RP_BNDRY_PLANES_CUT_NEWLINE}"
        "-DBNDRY_PLANES_EARLY_STOP=${ADD_TEST_RP_BNDRY_PLANES_EARLY_STOP}"
        "-DBNDRY_PLANES_READ_INPUT=${ADD_TEST_RP_BNDRY_PLANES_READ_INPUT}"
        "-DCHECKER=${ADD_TEST_RP_CHECKER}"
        "-DCHECK_MODE=${ADD_TEST_RP_CHECK_MODE}"
        "-DCHECKER_LOG=${_checker_log}"
        -P ${PROJECT_SOURCE_DIR}/Tests/RunRestartParity.cmake)
    # The reservation has to cover the widest leg, which need not be NP.
    set(_procs "${NP}")
    foreach(_n "${ADD_TEST_RP_CHK_NRANKS}" "${ADD_TEST_RP_RESTART_NRANKS}")
        if(NOT "${_n}" STREQUAL "" AND _n GREATER _procs)
            set(_procs "${_n}")
        endif()
    endforeach()

    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT ${_ctest_timeout}
        PROCESSORS ${_procs}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;restart-parity"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/straight/simulation.log;${CURRENT_TEST_BINARY_DIR}/restart/checkpoint.log;${CURRENT_TEST_BINARY_DIR}/restart/restart.log;${CURRENT_TEST_BINARY_DIR}/parity.log;${_checker_log};${_post_checkpoint_log}")
endfunction(add_test_restart_parity)

# FLARE's copy of the numeric log comparison (Tests/CompareDataLogs.cmake, with the column
# groups its conductor logs use): a comparator that accepts everything passes every test
# that uses it, so it needs its own test. Pure CMake, no run, hence the "unit" label.
add_test(CompareDataLogs_SelfTest ${CMAKE_COMMAND}
    "-DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/test_files/CompareDataLogs_SelfTest"
    -P ${PROJECT_SOURCE_DIR}/Tests/CompareDataLogsSelfTest.cmake)
set_tests_properties(CompareDataLogs_SelfTest PROPERTIES TIMEOUT 60 PROCESSORS 1 LABELS "unit")

#=============================================================================
# FLARE's fire helpers
#=============================================================================
# Fire smoke test: one deck of a fire suite under Exec/RegTests, run for NSTEPS
# steps on the regression rank count, passing when the run exits cleanly and
# the fire plotfile of the last step is written. The suites carry no gold
# files: their physics checks live in the run_*.sh scripts beside the decks,
# which run every variant to its stop time and are too long for CI. The whole
# suite directory is copied so a deck finds its inputs_base, sounding, fuel
# map, building list and schedules.
function(add_test_fire TEST_NAME SUITE_DIR INPUT_FILE NSTEPS)
    set(options )
    set(oneValueArgs "RUNTIME_OPTIONS" "NRANKS")
    set(multiValueArgs )
    cmake_parse_arguments(ADD_TEST_FIRE "${options}" "${oneValueArgs}"
        "${multiValueArgs}" ${ARGN})

    set(CURRENT_TEST_SOURCE_DIR ${PROJECT_SOURCE_DIR}/Exec/RegTests/${SUITE_DIR})
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_NAME})
    file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR})
    file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*")
    file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/")

    # NRANKS overrides the regression rank count: a 20-cell coarse deck has no
    # two-rank decomposition whose box edges divide by the fire grid ratio
    if(ERF_ENABLE_MPI)
        if("${ADD_TEST_FIRE_NRANKS}" STREQUAL "")
            set(NP ${ERF_TEST_NRANKS})
        else()
            set(NP ${ADD_TEST_FIRE_NRANKS})
        endif()
        set(MPI_COMMANDS "${ERF_SH_MPI_LAUNCHER}${MPIEXEC_NUMPROC_FLAG} ${NP} ${MPIEXEC_PREFLAGS}")
    else()
        set(NP 1)
        unset(MPI_COMMANDS)
    endif()

    resolve_test_exe("" "${ERF_TEST_EXE}" TEST_EXE)

    # fire plotfile names carry the step number padded to five digits
    set(_step "0000${NSTEPS}")
    string(LENGTH "${_step}" _len)
    math(EXPR _start "${_len} - 5")
    string(SUBSTRING "${_step}" ${_start} 5 _step)
    set(PLTFILE "plt_fire_${_step}")

    set(RUNTIME_OPTIONS "max_step=${NSTEPS} erf.fire.plot_int=${NSTEPS} erf.fire.plot_file=plt_fire_ erf.plot_int_1=-1 erf.plot_int_2=-1 erf.check_int=-1 ${ADD_TEST_FIRE_RUNTIME_OPTIONS}")
    set(test_log "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log")
    # the fire and dust CSVs are appended to (a restart keeps its earlier leg), so a
    # rerun in the same directory starts from a clean slate
    set(test_command sh -c "cd ${CURRENT_TEST_BINARY_DIR} && rm -rf plt_fire_* plt_dust_* fire_stats* dust_diag* && ${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${INPUT_FILE} ${RUNTIME_OPTIONS} > ${test_log} 2>&1 && test -f ${CURRENT_TEST_BINARY_DIR}/${PLTFILE}/Header")

    add_test(${TEST_NAME} ${test_command})
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 1800
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;fire"
        ATTACHED_FILES_ON_FAIL "${test_log}"
    )
endfunction(add_test_fire)

# Fire suite deck followed by its check script (pure Python; the plotfile reader
# erf_plotfile.py from Exec/Tools is copied next to it). The script's exit
# code is the verdict and its table is echoed into the ctest output.
function(add_test_fire_check TEST_NAME SUITE_DIR INPUT_FILE NSTEPS CHECK_SCRIPT)
    set(options )
    set(oneValueArgs "RUNTIME_OPTIONS" "NRANKS")
    set(multiValueArgs )
    cmake_parse_arguments(ADD_TEST_FIRE_CHECK "${options}" "${oneValueArgs}"
        "${multiValueArgs}" ${ARGN})

    set(CURRENT_TEST_SOURCE_DIR ${PROJECT_SOURCE_DIR}/Exec/RegTests/${SUITE_DIR})
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_NAME})
    file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR})
    file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*")
    file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/")
    file(COPY ${ERF_SOURCE_DIR}/Exec/Tools/erf_plotfile.py
         DESTINATION "${CURRENT_TEST_BINARY_DIR}/")

    if(ERF_ENABLE_MPI)
        if("${ADD_TEST_FIRE_CHECK_NRANKS}" STREQUAL "")
            set(NP ${ERF_TEST_NRANKS})
        else()
            set(NP ${ADD_TEST_FIRE_CHECK_NRANKS})
        endif()
        set(MPI_COMMANDS "${ERF_SH_MPI_LAUNCHER}${MPIEXEC_NUMPROC_FLAG} ${NP} ${MPIEXEC_PREFLAGS}")
    else()
        set(NP 1)
        unset(MPI_COMMANDS)
    endif()

    resolve_test_exe("" "${ERF_TEST_EXE}" TEST_EXE)
    set(RUNTIME_OPTIONS "max_step=${NSTEPS} erf.plot_int_1=-1 erf.plot_int_2=-1 erf.check_int=-1 ${ADD_TEST_FIRE_CHECK_RUNTIME_OPTIONS}")
    set(test_log "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log")
    set(check_log "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.check.log")
    # stale plotfiles and diagnostics of an earlier run are removed first: the
    # dust and fire CSVs are appended to (a restart keeps its earlier leg)
    set(test_command sh -c "cd ${CURRENT_TEST_BINARY_DIR} && rm -rf plt_fire_* plt_dust_* dust_diag* fire_stats* CHECK_FAILED && ${MPI_COMMANDS} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${INPUT_FILE} ${RUNTIME_OPTIONS} > ${test_log} 2>&1 || ( tail -n 60 ${test_log} && false ) && ( ${ERF_RANS_PYTHON} ${CURRENT_TEST_BINARY_DIR}/${CHECK_SCRIPT} > ${check_log} 2>&1 || touch ${CURRENT_TEST_BINARY_DIR}/CHECK_FAILED ) && cat ${check_log} && test ! -f ${CURRENT_TEST_BINARY_DIR}/CHECK_FAILED")

    add_test(${TEST_NAME} ${test_command})
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 1800
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;fire"
        ATTACHED_FILES_ON_FAIL "${test_log};${check_log}"
    )
endfunction(add_test_fire_check)

# Fire suite driven by its own run script: copy SUITE_DIR (and erf_plotfile.py from
# Canonical_RANS) and run SCRIPT with the executable, MPIRUN set for NRANKS ranks. The
# script's exit status is the verdict; it runs the variants and their checks itself.
function(add_test_fire_script TEST_NAME SUITE_DIR SCRIPT)
    set(options )
    set(oneValueArgs "NRANKS")
    set(multiValueArgs )
    cmake_parse_arguments(ADD_TEST_FIRE_SCRIPT "${options}" "${oneValueArgs}"
        "${multiValueArgs}" ${ARGN})

    set(CURRENT_TEST_SOURCE_DIR ${PROJECT_SOURCE_DIR}/Exec/RegTests/${SUITE_DIR})
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_NAME})
    file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR})
    file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*")
    file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/")
    file(COPY ${ERF_SOURCE_DIR}/Exec/Tools/erf_plotfile.py
         DESTINATION "${CURRENT_TEST_BINARY_DIR}/")

    if("${ADD_TEST_FIRE_SCRIPT_NRANKS}" STREQUAL "")
        set(NP ${ERF_TEST_NRANKS})
    else()
        set(NP ${ADD_TEST_FIRE_SCRIPT_NRANKS})
    endif()
    # unquoted on purpose: the scripts expand ${MPIRUN} unquoted, so ERF_SH_MPI_LAUNCHER's
    # quotes would reach the shell as literal characters; the callers are NOT WIN32 only
    set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NP} ${MPIEXEC_PREFLAGS}")

    resolve_test_exe("" "${ERF_TEST_EXE}" TEST_EXE)
    set(test_log "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log")
    # the script removes its own earlier output before running
    set(test_command sh -c "cd ${CURRENT_TEST_BINARY_DIR} && MPIRUN='${MPI_COMMANDS}' PYTHON=${ERF_RANS_PYTHON} FCOMPARE=${FCOMPARE_EXE} sh ${CURRENT_TEST_BINARY_DIR}/${SCRIPT} ${TEST_EXE} > ${test_log} 2>&1 || ( cat ${test_log} && false ) && cat ${test_log}")

    add_test(${TEST_NAME} ${test_command})
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 1800
        PROCESSORS ${NP}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;fire"
        ATTACHED_FILES_ON_FAIL "${test_log}"
    )
endfunction(add_test_fire_script)

# Fire deck that writes a checkpoint at step 1 on WRITE_NRANKS ranks and restarts
# from it on RESTART_NRANKS with RESTART_OPTIONS; passes only when the restart
# stops with EXPECTED_MESSAGE.
function(add_test_fire_restart_abort TEST_NAME SUITE_DIR INPUT_FILE WRITE_NRANKS RESTART_NRANKS EXPECTED_MESSAGE RESTART_OPTIONS)
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_NAME})
    file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR})
    file(GLOB TEST_FILES "${PROJECT_SOURCE_DIR}/Exec/RegTests/${SUITE_DIR}/*")
    file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/")
    set(MPI_WRITE "${ERF_SH_MPI_LAUNCHER}${MPIEXEC_NUMPROC_FLAG} ${WRITE_NRANKS} ${MPIEXEC_PREFLAGS}")
    set(MPI_RESTART "${ERF_SH_MPI_LAUNCHER}${MPIEXEC_NUMPROC_FLAG} ${RESTART_NRANKS} ${MPIEXEC_PREFLAGS}")
    resolve_test_exe("" "${ERF_TEST_EXE}" TEST_EXE)
    set(test_log "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log")
    set(write_log "${CURRENT_TEST_BINARY_DIR}/write.log")
    set(common "erf.plot_int_1=-1 erf.plot_int_2=-1 erf.fire.plot_int=-1 amrex.call_addr2line=0")
    set(test_command sh -c "rm -rf ${CURRENT_TEST_BINARY_DIR}/chk_abort* && ${MPI_WRITE} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${INPUT_FILE} max_step=1 erf.check_int=1 erf.check_file=chk_abort ${common} > ${write_log} 2>&1 && ${MPI_RESTART} ${TEST_EXE} ${CURRENT_TEST_BINARY_DIR}/${INPUT_FILE} max_step=2 erf.restart=chk_abort00001 erf.check_int=-1 ${common} ${RESTART_OPTIONS} 2>&1 | tee ${test_log}")
    add_test(${TEST_NAME} ${test_command})
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 600
        PROCESSORS ${RESTART_NRANKS}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;fire"
        PASS_REGULAR_EXPRESSION "${EXPECTED_MESSAGE}"
        ATTACHED_FILES_ON_FAIL "${test_log};${write_log}"
    )
endfunction(add_test_fire_restart_abort)

# Fire start-up check: add_test_abort on one deck of a fire suite under Exec/RegTests
function(add_test_fire_abort TEST_NAME SUITE_DIR INPUT_FILE EXPECTED_MESSAGE RUNTIME_OPTIONS)
    add_test_abort(${TEST_NAME} ${PROJECT_SOURCE_DIR}/Exec/RegTests/${SUITE_DIR} ${INPUT_FILE}
                   "${EXPECTED_MESSAGE}" "${RUNTIME_OPTIONS}")
    set_tests_properties(${TEST_NAME} PROPERTIES LABELS "regression;fire")
endfunction(add_test_fire_abort)

# Fire fuel map row order: run INPUT_FILE of a fire suite under Exec/RegTests to
# its step-0 fire plotfile and check fire_fuel_load at the south and north edges
# along x = X with amrex_fextract (Tests/RunFireFuelMapRows.cmake). The deck must
# load its ESRI ASCII map with erf.fire.fuel_map.load_from_map.
function(add_test_fire_fuel_map_rows TEST_NAME SUITE_DIR INPUT_FILE X SOUTH_MIN SOUTH_MAX NORTH_MIN NORTH_MAX)
    set(CURRENT_TEST_SOURCE_DIR ${PROJECT_SOURCE_DIR}/Exec/RegTests/${SUITE_DIR})
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_NAME})
    file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR})
    file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*")
    file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/")

    resolve_test_exe("" "${ERF_TEST_EXE}" TEST_EXE)
    string(REPLACE "amrex_fcompare" "amrex_fextract" FEXTRACT_EXE "${FCOMPARE_EXE}")

    add_test(${TEST_NAME} ${CMAKE_COMMAND}
        -DMPIEXEC=${MPIEXEC_EXECUTABLE}
        -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
        -DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}
        -DNRANKS=${ERF_TEST_NRANKS}
        -DTEST_EXE=${TEST_EXE}
        -DFEXTRACT=${FEXTRACT_EXE}
        -DINPUT=${CURRENT_TEST_BINARY_DIR}/${INPUT_FILE}
        -DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}
        -DX=${X}
        -DSOUTH_MIN=${SOUTH_MIN}
        -DSOUTH_MAX=${SOUTH_MAX}
        -DNORTH_MIN=${NORTH_MIN}
        -DNORTH_MAX=${NORTH_MAX}
        -P ${PROJECT_SOURCE_DIR}/Tests/RunFireFuelMapRows.cmake)
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 600
        PROCESSORS ${ERF_TEST_NRANKS}
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "regression;fire"
        ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/run.log")
endfunction(add_test_fire_fuel_map_rows)

# The fire plotfile table in fire_output.rst lists the fields "in this fixed
# order"; this checks it against the order ERF_FirePlotfileCatalog.H builds.
# Plain Python, no ERF run; the self-test proves a swapped or missing row fails.
if(ERF_TEST_PYTHON)
  add_test(NAME FirePlotfileDoc_MatchesCatalog
      COMMAND ${ERF_TEST_PYTHON} ${PROJECT_SOURCE_DIR}/Tests/check_fire_plotfile_doc.py
              --catalog ${PROJECT_SOURCE_DIR}/FireEngine/ERF_FirePlotfileCatalog.H
              --doc ${PROJECT_SOURCE_DIR}/Docs/sphinx_doc/theory/fire_output.rst)
  add_test(NAME FirePlotfileDoc_SelfTest
      COMMAND ${ERF_TEST_PYTHON} ${PROJECT_SOURCE_DIR}/Tests/check_fire_plotfile_doc.py --self-test)
  set_tests_properties(FirePlotfileDoc_MatchesCatalog FirePlotfileDoc_SelfTest
      PROPERTIES
      TIMEOUT 60
      PROCESSORS 1
      LABELS "unit;fire")
endif()

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    set(ERF_RANS_PYTHON "${Python3_EXECUTABLE}")
else()
    set(ERF_RANS_PYTHON "python3")
endif()

#=============================================================================
# Fire and dust smoke tests: one deck per suite under Exec/RegTests, a few
# steps each (ctest -L fire, or -R Fire)
#=============================================================================
if(ERF_ENABLE_FIRE)
add_test_fire(FireAccelerationClock_levelset_front FireAccelerationClock inputs_levelset_front 40)
add_test_fire(FireBurnout_base              FireBurnout           inputs_base                40)
add_test_fire(FireDirectionalShape_ellipse  FireDirectionalShape  inputs_ellipse             40)
add_test_fire(FireEmcModel_van_wagner       FireEmcModel          inputs_van_wagner          40)
add_test_fire(FireExposure_noib             FireExposure          inputs_noib                40)
add_test_fire(FireFbp_c2                    FireFbp               inputs_fbp_c2              40)
add_test_fire(FireFluxPartition_cfbm        FireFluxPartition     inputs_cfbm                40)
add_test_fire(FireHeatPlacement_add_noib    FireHeatPlacement     inputs_add_noib            40)
add_test_fire(FireHybridObstacles_noib      FireHybridObstacles   inputs_hybrid_noib         40)
add_test_fire(FireLevelSetEllipse_ellipse   FireLevelSetEllipse   inputs_ellipse             40)
add_test_fire(FireNearWall_noib_mask_wall   FireNearWall          inputs_noib_mask_wall      40)
add_test_fire(FirePerimeterIgnition_t0      FirePerimeterIgnition inputs_t0                  40)
add_test_fire(FireRestart_levelset_straight FireRestart           inputs_levelset_straight   40 NRANKS 1)
# the MRF fire thermal excess reads the lagged fire flux in the halo columns of every
# tile; the coupled restart deck with MRF and the option on stopped on an out-of-bound
# read in the first step before the flux carried ghost columns
add_test_fire(FireMrfThermalExcess          FireRestart           inputs_coupled_straight    40 NRANKS 1
    RUNTIME_OPTIONS "erf.pbl_type=MRF erf.pbl_mrf_fire_thermal_excess=true")
add_test_fire(FireRosComparison_rothermel   FireRosComparison     inputs_rothermel_isotropic 40 NRANKS 1)
# the defaults the 2026-10 validation changed, each next to the deck that keeps
# the earlier form: the grass model (Cheney, Gould and Catchpole 1998 against
# the fit it replaced), the BEHAVE net load (surface-area weighted against the
# summed), the FARSITE spread shape (Richards ellipse against the rectangle)
# and the wind limit (Rothermel's 0.9 I_R against the 300/500 ft/min rule)
add_test_fire(FireRosComparison_cheney_gould FireRosComparison    inputs_cheney_gould_isotropic 40 NRANKS 1)
add_test_fire(FireRosComparison_grass_simple FireRosComparison    inputs_grass_simple_isotropic 40 NRANKS 1)
add_test_fire(FireRosComparison_behave      FireRosComparison     inputs_behave_isotropic    40 NRANKS 1)
add_test_fire(FireRosComparison_behave_sum  FireRosComparison     inputs_behave_isotropic    40 NRANKS 1
    RUNTIME_OPTIONS "erf.fire.behave.net_load=sum")
add_test_fire(FireFarsiteShape_ellipse      FarsiteFrontUpdate    inputs_front_cell          40 NRANKS 1)
add_test_fire(FireFarsiteShape_rectangle    FarsiteFrontUpdate    inputs_rectangle           40 NRANKS 1)
add_test_fire(FireLineFire_wind5_cap        FireLineFire          inputs_wind5_cap           40 NRANKS 1)
add_test_fire(FireLineFire_wind5_cap_fuel_class FireLineFire      inputs_wind5_cap_fuel_class 40 NRANKS 1)
add_test_fire(FireScottBurgan_gr2           FireScottBurgan       inputs_sb_gr2              40)
add_test_fire(FireCustomFuel_uniform        FireCustomFuel        inputs_custom_grass        40)
add_test_fire(FireCustomFuel_map            FireCustomFuel        inputs_custom_map          40)
# erf.fire.rothermel_per_fuel with the level-set path: the deck differs from
# FireCustomFuel_map only in erf.fire.fuel_model_id, which a per-fuel run over a
# load_from_map raster must ignore entirely -- the Rothermel coefficients, the
# wind adjustment factor built from the bed depth and the Byram diagnostics built
# from the initial load each used to read it. run_custom_fuel.sh compares the two
# runs line for line and, with fcompare, field for field.
add_test_fire(FireCustomFuel_map_altid      FireCustomFuel        inputs_custom_map_altid    40)
# a custom code the deck declares non-burnable needs no properties block
add_test_fire(FireCustomFuel_nonburnable    FireCustomFuel        inputs_undeclared_nonburnable 40)
add_test_fire(FireStickMoisture_stick       FireStickMoisture     inputs_stick               40)
add_test_fire(FireStructureIgnition_on      FireStructureIgnition inputs_on                  40 NRANKS 1)
add_test_fire(FirePrecipSource_atmosphere   FirePrecipSource      inputs_atmosphere          40 NRANKS 1)
add_test_fire(FireLiveMoisture_fixed        FireLiveMoisture      inputs_fixed               40)
add_test_fire(FireWindSampling_sample20     FireWindSampling      inputs_sample20            40)
add_test_fire(FirePrescribed_ros_circle     FirePrescribed        inputs_ros_circle          40)
add_test_fire(FirePrescribed_heat_patch     FirePrescribed        inputs_heat_patch          40)
add_test_fire(FireFarsiteDefault            FarsiteDefault        inputs                     40)
add_test_fire(FireLevelSetPropagation       LevelSetPropagation   inputs                     40)
# fire without a surface layer at zlo must stop at start-up, not crash in the first step
add_test_fire_abort(FireNoSurfaceLayer_abort  FireRestart           inputs_levelset_straight
    "requires a surface layer at the bottom boundary" "zlo.type=SlipWall")
# a deck-defined fuel model has to be described completely and in range: a code
# outside 1000-1015, a missing property, a bed depth under the floor where the
# Balbi models silently return zero spread, a heat content left in BTU/lb, and a
# raster code no block defines all have to stop the run at start-up
add_test_fire_abort(FireCustomFuelBadCode_abort     FireCustomFuel inputs_bad_code
    "outside the custom range 1000-1015" "max_step=1")
add_test_fire_abort(FireCustomFuelMissing_abort     FireCustomFuel inputs_bad_missing
    "erf.fire.custom_fuel.1000.sav_1h_1_m is required" "max_step=1")
add_test_fire_abort(FireCustomFuelBadDepth_abort    FireCustomFuel inputs_bad_depth
    "erf.fire.custom_fuel.1000.depth_m must be in" "max_step=1")
add_test_fire_abort(FireCustomFuelBadHeat_abort     FireCustomFuel inputs_bad_heat
    "erf.fire.custom_fuel.1000.heat_content_J_kg must be in" "max_step=1")
# a 1-h SAV below 100 1/ft aborts; Rothermel's reaction term clamped it up to
# that, silently, until November 2026
add_test_fire_abort(FireCustomFuelLowSav_abort      FireCustomFuel inputs_custom_heavy
    "erf.fire.custom_fuel.1000.sav_1h_1_m must be in .328.084" "erf.fire.custom_fuel.1000.sav_1h_1_m=200.0 max_step=1")
add_test_fire_abort(FireCustomFuelBurnout_abort     FireCustomFuel inputs_bad_burnout
    "needs erf.fire.custom_fuel.1000.burnout_time_s" "max_step=1")
add_test_fire_abort(FireCustomFuelUndeclared_abort  FireCustomFuel inputs_bad_undeclared
    "the fuel map holds code 1007" "max_step=1")
add_test_fire_abort(FireCustomFuelUniformCode_abort FireCustomFuel inputs_bad_uniform_code
    "erf.fire.fuel_model_id = 1000 is in the custom range" "max_step=1")
# a fuel map is placed by cell index, so it must have the fire grid's size
add_test_fire_abort(FireFuelMapSize_abort     FireScottBurgan       inputs_sb_map
    "has 256 x 128 cells but the fire grid has 128 x 64" "erf.fire.grid_ratio=2")
# the first data row of an ESRI ASCII fuel map is the north edge: the sb_map deck's
# map puts TL3 (1.2329 kg/m2) along the north and NB8 water (no fuel) along the south
# execute_process needs mpiexec, and does not expand the executable globs used on Windows
if(ERF_ENABLE_MPI AND NOT WIN32)
add_test_fire_fuel_map_rows(FireScottBurgan_map_rows FireScottBurgan inputs_sb_map 100
    -0.000001 0.000001 1.2329 1.2330)
endif()
# a misspelt selector or a value the kernels cannot use stops at start-up
add_test_fire_abort(FireBadRosModel_abort     FireRestart           inputs_levelset_straight
    "erf.fire.ros_model = \"rothermal\" is not one of" "erf.fire.ros_model=rothermal")
add_test_fire_abort(FireBadCoupling_abort     FireRestart           inputs_levelset_straight
    "erf.fire.coupling_type = \"laged\" is not one of" "erf.fire.coupling_type=laged")
# selectors and ranges of the keys the 2026-10 validation added
add_test_fire_abort(FireBadWindLimit_abort    FireRestart           inputs_levelset_straight
    "erf.fire.wind_limit = \"rothermal\" is not one of" "erf.fire.wind_limit=rothermal")
add_test_fire_abort(FireBadFarsiteShape_abort FireRestart           inputs_levelset_straight
    "erf.fire.farsite.shape = \"oval\" is not one of" "erf.fire.farsite.shape=oval")
add_test_fire_abort(FireBadNetLoad_abort      FireRestart           inputs_levelset_straight
    "erf.fire.behave.net_load = \"average\" is not one of" "erf.fire.ros_model=behave erf.fire.behave.net_load=average")
add_test_fire_abort(FireCheneyKeyNoModel_abort FireRestart          inputs_levelset_straight
    "erf.fire.cheney_gould.pasture needs erf.fire.ros_model = cheney_gould" "erf.fire.cheney_gould.pasture=grazed")
# unit guards (November 2026 unit audit): a value in the wrong unit aborts with the
# input's name instead of running (foliar moisture in percent, temperatures in
# degrees C), and the two removed options say what replaces them
add_test_fire_abort(FireFoliarMoisturePercent_abort FireRestart     inputs_levelset_straight
    "erf.fire.crown.foliar_moisture is a fraction .1.0 = 100 %. and must be in .0.3, 3.0., got 100" "erf.fire.crown.foliar_moisture=100")
add_test_fire_abort(FireFlameTempCelsius_abort FireRestart          inputs_levelset_straight
    "erf.fire.flame_temp_T_amb is in K and must be in .150, 400., got 25" "erf.fire.flame_temp_T_amb=25")
add_test_fire_abort(FireFlameTiltCelsius_abort FireRestart          inputs_levelset_straight
    "erf.fire.flame_tilt_T_amb is in K and must be in .150, 400., got 25" "erf.fire.flame_tilt_T_amb=25")
add_test_fire_abort(FireThresholdTempCelsius_abort FireRestart      inputs_levelset_straight
    "erf.fire.ignition.threshold_temp is in K and must be in .305, 2000., got 300" "erf.fire.ignition.threshold_temp=300")
add_test_fire_abort(FireWafBehaviorplus_abort  FireRestart          inputs_levelset_straight
    "erf.fire.waf_formula = .behaviorplus. was removed in November 2026" "erf.fire.waf_formula=behaviorplus")
add_test_fire_abort(FireFlameTempNelson_abort  FireRestart          inputs_levelset_straight
    "erf.fire.flame_temp_method = .nelson_emissivity. was removed in November 2026" "erf.fire.flame_temp_method=nelson_emissivity")
add_test_fire_abort(FireWindZ0Above20ft_abort  FireRestart          inputs_levelset_straight
    "erf.fire.wind_sample_z0 must be below 6.096 m" "erf.fire.wind_ref_ht=20 erf.fire.wind_sample_z0=7")
add_test_fire_abort(FireWindZ0Profile_abort FireRestart inputs_levelset_straight
    "erf.fire.wind_sample_z0 must be at most 2.24" "erf.fire.wind_sample_z0=3")
add_test_fire_abort(FireSpotThermalModel_abort FireRestart          inputs_spotting_straight
    "erf.fire.spotting.thermal_fuel_model must be 0 or an Anderson model in Albini .1983. Table 4" "erf.fire.spotting.thermal_fuel_model=8")
add_test_fire_abort(FireSpotNoThermalFit_abort FireRestart          inputs_spotting_straight
    "has no Albini .1983. Table 4 thermal fit" "erf.fire.fuel_model_id=8")
add_test_fire_abort(FireSpotSb40NoThermal_abort FireRestart         inputs_spotting_straight
    "erf.fire.spotting.fuel_system = .40. needs erf.fire.spotting.thermal_fuel_model" "erf.fire.fuel_map.fuel_set=scott_burgan40 erf.fire.fuel_model_id=102")
add_test_fire_abort(FireSpotStructureNoThermal_abort FireStructureIgnition inputs_on
    "erf.fire.structures.ignition.enable with spotting needs erf.fire.spotting.thermal_fuel_model" "erf.fire.spotting.thermal_fuel_model=0")
add_test_fire_abort(FireSpotMapNoThermalFit_abort FireExposure  inputs_noib_spotting
    "erf.fire.fuel_map.file has fuel code 0 on cells with fuel" "erf.fire.fuel_map.nonburnable_codes=99")
add_test_fire_abort(FireFbpAndersonLw_abort    FireRestart          inputs_levelset_straight
    "erf.fire.directional_ellipse_lw = .anderson. takes the midflame wind" "erf.fire.ros_model=fbp erf.fire.directional_ros=true erf.fire.directional_shape=ellipse erf.fire.directional_ellipse_lw=anderson")
add_test_fire_abort(FireMacArthurCap_abort    FireRestart           inputs_levelset_straight
    "erf.fire.macarthur.ros_max must be >= 0" "erf.fire.macarthur.ros_max=-1")
# the model checks of the keys read only by one model
add_test_fire_abort(FireMacArthurCapNoModel_abort FireRestart       inputs_levelset_straight
    "erf.fire.macarthur.ros_max needs erf.fire.ros_model = macarthur" "erf.fire.ros_model=rothermel erf.fire.macarthur.ros_max=3")
add_test_fire_abort(FireNetLoadNoBehave_abort FireRestart           inputs_levelset_straight
    "erf.fire.behave.net_load needs erf.fire.ros_model = behave" "erf.fire.ros_model=rothermel erf.fire.behave.net_load=sum")
add_test_fire_abort(FireBadWindBelow_abort    FireRestart           inputs_levelset_straight
    "erf.fire.wind_below_first_cell = \"linear\" is not one of" "erf.fire.wind_below_first_cell=linear")
add_test_fire_abort(FireHybridGrassPair_abort FireRestart           inputs_levelset_straight
    "cheney_gould and grass_simple share one grass state" "erf.fire.ros_model=hybrid erf.fire.hybrid.primary=cheney_gould erf.fire.hybrid.secondary=grass_simple")
# a threshold ignition with no disc (erf.fire.ignition_r = 0) is a deferred ignition: the
# start-up must not abort for "no cells were marked as burned" (it did until 2026-10)
add_test_fire(FireThresholdNoDisc FireThresholdIgnition inputs_on 5 RUNTIME_OPTIONS "erf.fire.ignition_r=0")
# start-up checks added by the October 2026 audit: a value the kernels cannot use,
# or an input given without the switch that reads it, stops the run naming the key
# (each of these ran on silently before: Balbi returned 15 m/s with no bisection
# step, code 50 burned as fuel model 1, tau_wind = 0 divided by zero, and the
# stick/suppression/thermal-excess keys were ignored)
add_test_fire_abort(FireBalbiMaxIter_abort    FireRestart           inputs_levelset_straight
    "erf.fire.balbi.max_iter must be >= 1" "erf.fire.balbi.max_iter=0")
add_test_fire_abort(FireFuelModelId_abort     FireRestart           inputs_levelset_straight
    "is not an Anderson code" "erf.fire.fuel_model_id=50")
add_test_fire_abort(FireAccelTauWind_abort    FireRestart           inputs_levelset_straight
    "erf.fire.accel.tau_wind must be > 0" "erf.fire.accel.enable=true erf.fire.accel.enable_wind_lag=true erf.fire.accel.tau_wind=0")
add_test_fire_abort(FireStickKeyNoStick_abort FireRestart           inputs_levelset_straight
    "need erf.fire.moisture_model = stick" "erf.fire.moisture_dynamic=true erf.fire.stick.n_shells=4")
add_test_fire_abort(FireSuppressionKeyNoEnable_abort FireRestart    inputs_levelset_straight
    "needs erf.fire.suppression.enable = true" "erf.fire.suppression.file=actions.txt")
# the default source_mode = overwrite would discard another theta source on the fire
# level (here the Rayleigh damping); the deck must choose add or overwrite explicitly
add_test_fire_abort(FireSourceModeConflict_abort FireRestart        inputs_coupled_straight
    "the default erf.fire.source_mode = overwrite replaces the theta source" "erf.rayleigh_damp_T=true")
add_test_fire_abort(FireMrfExcessNoFlux_abort FireRestart           inputs_levelset_straight
    "erf.pbl_mrf_fire_thermal_excess = true on level 0 needs erf.fire.enable" "erf.pbl_type=MRF erf.pbl_mrf_fire_thermal_excess=true")
# the fire engine behind the shared coupling layer (erf.fire.engine): a name that is
# not an engine stops the run naming the key and the engines there are
add_test_fire_abort(FireEngineUnknown_abort   FireRestart           inputs_levelset_straight
    "erf.fire.engine = .lsfire. is not one of .*.flare." "erf.fire.engine=lsfire")
# 20 cells with max_grid_size 8 and blocking factor 1 chop into 7 + 7 + 6, none a
# multiple of the grid ratio: the run stops naming the box. For the power-of-2 ratio
# 4 the advice is a blocking factor; for 5 (the default ratio), which no AMReX
# blocking factor can be a multiple of, it is a max_grid_size
add_test_fire_abort(FireBoxRatio_abort        FireRestart           inputs_levelset_straight
    "Box sizes not divisible by grid_ratio: box 0 of level 0 is 7 x 7 cells.*Set amr.blocking_factor_x"
    "amr.max_grid_size_x=8 amr.max_grid_size_y=8")
add_test_fire_abort(FireBoxRatio5_abort       FireRestart           inputs_levelset_straight
    "Box sizes not divisible by grid_ratio: box 0 of level 0 is 7 x 7 cells.*which 5 is not, so set the blocking factor to 1"
    "erf.fire.grid_ratio=5 amr.max_grid_size_x=8 amr.max_grid_size_y=8")
# a restart that re-makes the level-0 grids stops before the fire grid is built or
# any checkpoint field is read; the fire used to trip its grid-ratio check or VisMF's
# BoxArray check part-way through. Both triggers: erf.regrid_level_0_on_restart, and
# a restart on more ranks than the checkpoint has level-0 boxes (one box, two ranks)
if(ERF_ENABLE_MPI AND NOT WIN32)
add_test_fire_restart_abort(FireRestartRegridLevel0_abort FireRestart inputs_levelset_straight 1 1
    "the fire model .its fire grid is built from the run's level boxes"
    "erf.regrid_level_0_on_restart=1 amr.max_grid_size_x=8 amr.max_grid_size_y=8")
add_test_fire_restart_abort(FireRestartMoreRanks_abort FireRestart inputs_levelset_straight 1 2
    "the fire model .its fire grid is built from the run's level boxes" "")
# the lagged fire flux the coupling layer saves in the checkpoint: a restarted coupled
# run gives the atmosphere and fire plotfiles of a straight run at zero tolerance, on
# nine boxes over two ranks (a restart that dropped the flux injects no fire heat in
# its first step)
add_test_fire_script(FireRestart_coupled     FireRestart           run_coupled_restart.sh NRANKS 2)
endif()
# the fire grid on a refined level (erf.fire.anchor_level, the finest level by default):
# the same front as a single-level run at that resolution, the level-0 heat budget after
# average-down (and its loss with the fire on the coarser level), restart, and a restart
# that would move the fire grid stopping at start-up
if(ERF_ENABLE_MPI AND NOT WIN32)
add_test_fire_script(FireAnchorLevel          FireAnchorLevel       run_anchor_level.sh NRANKS 2)
endif()
# where the rain that wets the dead fuel comes from (erf.fire.precip_source): the
# atmosphere's rain per column wets the fuel under the raining columns only and the
# fire-grid rate equals the change of the Kessler surface accumulation over the step,
# the uniform rate wets every cell, and the restart carries the accumulation snapshot
if(ERF_ENABLE_MPI AND NOT WIN32)
add_test_fire_script(FirePrecipSource         FirePrecipSource      run_precip.sh NRANKS 1)
endif()
# its start-up checks: a scheme without precipitation accumulators, static moisture,
# and both rain sources named at once
add_test_fire_abort(FirePrecipSource_norain_abort  FirePrecipSource inputs_atmosphere
    "provides none" "erf.moisture_model=Kessler_NoRain")
add_test_fire_abort(FirePrecipSource_static_abort  FirePrecipSource inputs_atmosphere
    "holds them fixed" "erf.fire.moisture_dynamic=false")
add_test_fire_abort(FirePrecipSource_two_sources_abort FirePrecipSource inputs_atmosphere
    "the rain has one source" "erf.fire.precip_rate_mm_hr=1.0")
# its start-up checks: a level above the finest, a regridding level, a refinement box
# short of the domain top, two separate patches, and the dust layer (level 0 only)
add_test_fire_abort(FireAnchorLevel_above_finest_abort FireAnchorLevel inputs_base
    "is above the finest level of this run" "erf.fire.anchor_level=2")
add_test_fire_abort(FireAnchorLevel_regrid_abort       FireAnchorLevel inputs_base
    "regrids it" "erf.regrid_int=10")
add_test_fire_abort(FireAnchorLevel_partial_height_abort FireAnchorLevel inputs_partial_height
    "Cannot decompose in z direction" "")
add_test_fire_abort(FireAnchorLevel_two_patches_abort  FireAnchorLevel inputs_two_patches
    "but the fire grid needs one rectangle" "")
if(ERF_ENABLE_DUST)
add_test_fire_abort(FireAnchorLevel_dust_abort         FireAnchorLevel inputs_base
    "The dust layer and the fire-dust coupling run on level 0" "erf.dust.enable=true")
endif()
# the fire at the edge of the fire grid: the guard band records the first contact in the
# statistics CSV and warns (warn), never fires on a fire far from every wall (far), or
# stops the run on a disc that starts inside the band (abort)
# The three ways fires meet (FireMergingFronts): spot fires coalescing from the
# schedule, two lines meeting at an angle (a junction fire, one V polyline) and two
# parallel lines from two perimeter files (erf.fire.ignition.polygon_file takes a
# list). Still air, a prescribed 1 m/s, so every arrival time follows from the
# distance to the nearest ignition; the checker compares the whole field and the
# cells where the fronts meet, and that every ignition is stamped (with one file
# read, as the code did before the list, the north line of the parallel pair is
# unburned and the strip between the lines arrives 20 s late).
add_test_fire_check(FireMergingFronts_coalescing FireMergingFronts inputs_coalescing 40 check_merging_fronts.py)
add_test_fire_check(FireMergingFronts_junction   FireMergingFronts inputs_junction   40 check_merging_fronts.py)
add_test_fire_check(FireMergingFronts_parallel   FireMergingFronts inputs_parallel   40 check_merging_fronts.py)
# erf.fire.farsite_winds: FARSITE/FlamMap gridded winds (_vel/_ang ESRI ASCII
# grids and a .atm file, the layout WindNinja writes), on two ranks over three
# boxes so the gather sees an uneven split. Prescribed u = 3, v = 4 m/s gives
# every cell in closed form (11.18 mph from 217 deg); the atmospheric deck's
# 2 x 2 block means are checked against fire_wind_ref in the step-61 plotfile.
add_test_fire_check(FireFarsiteWinds_prescribed FireFarsiteWinds inputs_prescribed 125 check_farsite_winds.py NRANKS 2)
add_test_fire_check(FireFarsiteWinds_atmosphere FireFarsiteWinds inputs_atmosphere 125 check_farsite_winds.py NRANKS 2)
# a restart from the step-61 checkpoint drops the straight run's later .atm
# line before appending, so the file and the grids match the straight run
if(NOT WIN32)
add_test_fire_script(FireFarsiteWinds_restart FireFarsiteWinds run_restart.sh NRANKS 2)
endif()
# the .atm units fix the wind height (ENGLISH 20 ft, METRIC 10 m), the dates
# need start_datetime, the interval is in whole minutes, and the keys need it on
add_test_fire_abort(FireFarsiteWinds_height_abort FireFarsiteWinds inputs_prescribed
    "erf.fire.farsite_winds.units = metric means a wind at 10.000000 m" "erf.fire.farsite_winds.units=metric")
add_test_fire_abort(FireFarsiteWinds_no_start_abort FireFarsiteWinds inputs_no_start
    "erf.fire.farsite_winds.interval > 0 needs start_datetime" "")
add_test_fire_abort(FireFarsiteWinds_interval_abort FireFarsiteWinds inputs_prescribed
    "erf.fire.farsite_winds.interval must be 0 .off. or >= 60 s" "erf.fire.farsite_winds.interval=30")
add_test_fire_abort(FireFarsiteWinds_key_without_interval_abort FireFarsiteWinds inputs_key_only
    "erf.fire.farsite_winds.units needs erf.fire.farsite_winds.interval > 0" "")
add_test_fire_abort(FireFarsiteWinds_coarsen_abort FireFarsiteWinds inputs_atmosphere
    "erf.fire.farsite_winds.coarsen = 7 must divide the fire grid, 60 x 20 cells" "erf.fire.farsite_winds.coarsen=7")
# every file of the list is read: a second file that does not exist stops the run
add_test_fire_abort(FireMergingFronts_missing_file_abort FireMergingFronts inputs_missing_file
    "Cannot open polygon vertex file: parallel_missing.csv" "")
# a polyline file with one vertex has no segment and would mark nothing: stop naming it
add_test_fire_abort(FireMergingFronts_one_vertex_abort FireMergingFronts inputs_one_vertex
    "'one_vertex.csv' has 1 vertices; a polyline needs at least 2" "")
# erf.fire.wind_sample_z0 must stay below the height of the reference wind:
# wind_ref_ht, or 6.096 m with use_per_fuel_wind_ht (until 2026-10 the per-fuel
# case was checked against wind_ref_ht, so z0 = 8 with wind_ref_ht = 10 ran with
# the ember drift's log profile reversed)
add_test_fire_abort(FireWindSampleZ0_abort FireLineFire inputs_nowind
    "erf.fire.wind_sample_z0 must be positive and below the height of the reference wind" "erf.fire.wind_sample_z0=7")
add_test_fire_abort(FireWindSampleZ0PerFuel_abort FireLineFire inputs_nowind
    "erf.fire.wind_sample_z0 must be positive and below the height of the reference wind.*with that height 6.096" "erf.fire.wind_ref_ht=10 erf.fire.wind_sample_z0=8 erf.fire.use_per_fuel_wind_ht=true")
# use_per_fuel_wind_ht samples the wind at 6.096 m without a fuel map too (until
# 2026-10 a uniform-fuel run sampled at wind_ref_ht while the checks and the
# ember drift took 6.096 m; this deck aborted, 8 m not being above 10 m)
add_test_fire_abort(FireWindPerFuelUniform_height FireLineFire inputs_wind2p5
    "Wind sampled at 8 m above ground, log-law factor to 6.096" "erf.fire.wind_ref_ht=10 erf.fire.wind_sample_ht=8 erf.fire.use_per_fuel_wind_ht=true max_step=1")
# the same without resampling: the height the extraction records in every
# fire cell (the old code sampled at 10 m), and the warning that wind_ref_ht is unused
add_test_fire_abort(FireWindPerFuelUniform_extract FireLineFire inputs_wind2p5
    "Wind extraction height range: min=6.096 m  max=6.096 m" "erf.fire.wind_ref_ht=10 erf.fire.use_per_fuel_wind_ht=true max_step=1")
add_test_fire_abort(FireWindRefHtPerFuel_warn FireLineFire inputs_wind2p5
    "WARNING: erf.fire.wind_ref_ht = 10 m is not used" "erf.fire.wind_ref_ht=10 erf.fire.use_per_fuel_wind_ht=true max_step=1")
# the start-up summary of the models in use is printed with fire_debug off, and
# marks a key the inputs set and one they do not
add_test_fire_abort(FireModelSummary_print FireLineFire inputs_wind2p5
    "Models \\(erf.fire key = value.*reaction_velocity_formula = albini  \\(default\\).*wind limit form +wind_limit = fuel_class  \\(set\\)" "erf.fire.fire_debug=false erf.fire.use_wind_limit=true erf.fire.wind_limit=fuel_class max_step=0")
# the domain must be deeper than the wind's sampling height; until 2026-10 the
# check compared the absolute top, so a 4 m deep domain at a raised bottom (here
# z = 145 to 149 m) ran.
# Since erf-model/ERF#4827 the sounding starts at geometry.prob_lo[2], hence the
# _z145 soundings (the same profiles shifted up to the raised bottom)
add_test_fire_abort(FireDomainDepth_abort FireLineFire inputs_nowind
    "Domain depth 4.0+ m .* should exceed the wind's sampling height 6.1" "geometry.prob_lo=0\ 0\ 145 geometry.prob_hi=400\ 80\ 149 erf.input_sounding_file=input_sounding_u0_z145 max_step=1")
# with erf.fire.wind_sample_ht the wind is sampled there: a 15 m deep domain
# cannot sample at 20 m (the check compared the 6.1 m target until 2026-10)
add_test_fire_abort(FireDomainDepthSample_abort FireLineFire inputs_wind2p5
    "Domain depth 15.0+ m .* should exceed the wind's sampling height 20" "geometry.prob_lo=0\ 0\ 145 geometry.prob_hi=400\ 80\ 160 erf.input_sounding_file=input_sounding_u2p5_z145 erf.fire.wind_sample_ht=20 max_step=1")

# A firebreak with a t = 0 polyline ignition stays in the non-burnable mask: all
# 750 break cells. Until 2026-10 the mask read the break back from phi after the
# polyline's whole-grid distance had been merged over it, and came out empty.
add_test_fire_abort(FireMergingFronts_firebreak_mask FireMergingFronts inputs_firebreak
    "Non-burnable mask: 750 cells" "erf.fire.fire_debug=1")
# Suppression (erf.fire.suppression.*): each scenario on the level-set and the
# FARSITE path, checked from the last fire plotfile and the suppression log
add_test_fire_check(FireSuppression_line_early_levelset FireSuppression inputs_line_early 40 check_suppression.py NRANKS 1)
add_test_fire_check(FireSuppression_line_early_farsite  FireSuppression inputs_line_early 40 check_suppression.py NRANKS 1
                    RUNTIME_OPTIONS "erf.fire.propagation_method=farsite")
add_test_fire_check(FireSuppression_line_late_levelset  FireSuppression inputs_line_late  40 check_suppression.py NRANKS 1)
add_test_fire_check(FireSuppression_line_late_farsite   FireSuppression inputs_line_late  40 check_suppression.py NRANKS 1
                    RUNTIME_OPTIONS "erf.fire.propagation_method=farsite")
add_test_fire_check(FireSuppression_drop_hold_levelset  FireSuppression inputs_drop_hold  40 check_suppression.py NRANKS 1)
add_test_fire_check(FireSuppression_drop_hold_farsite   FireSuppression inputs_drop_hold  40 check_suppression.py NRANKS 1
                    RUNTIME_OPTIONS "erf.fire.propagation_method=farsite")
add_test_fire_check(FireSuppression_drop_slow_levelset  FireSuppression inputs_drop_slow  40 check_suppression.py NRANKS 1)
add_test_fire_check(FireSuppression_drop_slow_farsite   FireSuppression inputs_drop_slow  40 check_suppression.py NRANKS 1
                    RUNTIME_OPTIONS "erf.fire.propagation_method=farsite")
add_test_fire_check(FireSuppression_hold_levelset       FireSuppression inputs_hold       40 check_suppression.py NRANKS 1)
add_test_fire_check(FireSuppression_hold_farsite        FireSuppression inputs_hold       40 check_suppression.py NRANKS 1
                    RUNTIME_OPTIONS "erf.fire.propagation_method=farsite")
add_test_fire_check(FireSuppression_burnout_levelset    FireSuppression inputs_burnout    40 check_suppression.py NRANKS 1)
add_test_fire_check(FireSuppression_burnout_farsite     FireSuppression inputs_burnout    40 check_suppression.py NRANKS 1
                    RUNTIME_OPTIONS "erf.fire.propagation_method=farsite")
add_test_fire_abort(FireSuppression_bad_line_abort      FireSuppression inputs_line_early
                    "suppression file 'actions_bad.txt' line 2: rate must be a number > 0"
                    "erf.fire.suppression.file=actions_bad.txt")
add_test_fire_abort(FireSuppression_duplicate_id_abort  FireSuppression inputs_line_early
                    "suppression file 'actions_duplicate.txt' line 3: duplicate id"
                    "erf.fire.suppression.file=actions_duplicate.txt")
add_test_fire_abort(FireSuppression_no_file_abort       FireSuppression inputs_base
                    "erf.fire.suppression.enable needs erf.fire.suppression.file"
                    "erf.v=0")
# The watched file (a writer process appends the line while the run polls, on
# one rank and on two, where the re-read is broadcast) and the restart with a
# line under construction, both on both propagation paths
if(ERF_ENABLE_MPI AND NOT WIN32)
add_test_fire_script(FireSuppression_poll_levelset    FireSuppression run_poll.sh NRANKS 1)
add_test_fire_script(FireSuppression_poll_farsite     FireSuppression run_poll_farsite.sh NRANKS 1)
add_test_fire_script(FireSuppression_poll_2ranks      FireSuppression run_poll.sh NRANKS 2)
add_test_fire_script(FireSuppression_restart_levelset FireSuppression run_restart.sh NRANKS 1)
add_test_fire_script(FireSuppression_restart_farsite  FireSuppression run_restart_farsite.sh NRANKS 1)
endif()

add_test_fire_check(FireBoundaryGuard_warn    FireBoundaryGuard     inputs_warn  40 check_guard.py NRANKS 1)
add_test_fire_check(FireBoundaryGuard_far     FireBoundaryGuard     inputs_far   40 check_guard.py NRANKS 1)
add_test_fire_abort(FireBoundaryGuard_abort   FireBoundaryGuard     inputs_abort
    "reached the boundary guard band" "")
# structure ignition and house-to-house spread (erf.fire.structures.ignition.*): the
# old and the new path on the same three houses, a radiation-only variant in which
# the second house can only ignite from the first, a restart, and the checker on the
# old path, where it must fail
if(ERF_ENABLE_MPI AND NOT WIN32)
add_test_fire_script(FireStructureIgnition   FireStructureIgnition run_structure_ignition.sh NRANKS 1)
endif()
# its start-up checks: ignition without the exposure accumulators it reads, and a
# burn curve whose growth phase alone would release more than the curve allows
add_test_fire_abort(FireStructureIgnition_no_exposure_abort FireStructureIgnition inputs_on
    "needs erf.fire.exposure.enable" "erf.fire.exposure.enable=false")
add_test_fire_abort(FireStructureIgnition_curve_abort FireStructureIgnition inputs_on
    "exceeds .1 - decay_fraction. of fuel_load_J_m2" "erf.fire.structures.ignition.growth_time_s=100.0")
add_test_fire_abort(FireStructureIgnition_decay_abort FireStructureIgnition inputs_on
    "decay_fraction must be in .0, 1.." "erf.fire.structures.ignition.decay_fraction=1.0")
# the Rothermel reaction-velocity exponent (erf.fire.reaction_velocity_formula) and
# the WRF-Fire fuel-load deflation (erf.fire.wrf_bmst_compat): the head rate of each
# combination, calm and windy, against its closed form, on the uniform-fuel path and
# through a fuel map's per-fuel table; and the checker's own pass/fail logic, where a
# run that ignored either option must fail
if(ERF_ENABLE_MPI AND NOT WIN32)
add_test_fire_script(FireReactionVelocityFormula FireReactionVelocityFormula run_regtest.sh NRANKS 1)
# the wind brought to each model's height (FBP: 10 m, the WAF: 20 ft) from the
# sampled height with the neutral log profile; FBP and the WAF read the sampled
# wind as their own height's until November 2026
add_test_fire_script(FireWindHeights          FireFbp               run_wind_heights.sh NRANKS 1)
endif()
add_test(FireReactionVelocityFormula_SelfTest ${ERF_RANS_PYTHON}
    ${PROJECT_SOURCE_DIR}/Exec/RegTests/FireReactionVelocityFormula/check_regtest.py --self-test)
set_tests_properties(FireReactionVelocityFormula_SelfTest
    PROPERTIES
    TIMEOUT 60
    PROCESSORS 1
    LABELS "fire;unit")
# both keys are ignored by a model that is not Rothermel, and the run says so,
# once: FireLayer's copy of the fire inputs used to read them again, so every
# start-up warning of FireParams was printed twice
add_test_fire_abort(FireReactionVelocityFormula_behave_warning FireReactionVelocityFormula inputs_rothermel_bmst
    "erf.fire.wrf_bmst_compat is ignored with erf.fire.ros_model = behave" "erf.fire.ros_model=behave")
set_tests_properties(FireReactionVelocityFormula_behave_warning PROPERTIES
    FAIL_REGULAR_EXPRESSION "wrf_bmst_compat is ignored.*wrf_bmst_compat is ignored")
# and a misspelt formula stops at start-up
add_test_fire_abort(FireReactionVelocityFormula_bad_abort FireReactionVelocityFormula inputs_albini
    "erf.fire.reaction_velocity_formula = \"rothermal\" is not one of" "erf.fire.reaction_velocity_formula=rothermal")
# level-set start-up checks (Exec/CanonicalTests/Fire/Unit_Tests/LevelSetStartupChecks):
# a misspelt reinit_clamp, and an artificial viscosity or cfl past the substep's
# linear stability limit for any front direction, stop the run
function(add_test_fire_levelset_startup TEST_NAME INPUT_FILE EXPECTED_MESSAGE)
    add_test_abort(${TEST_NAME} ${PROJECT_SOURCE_DIR}/Exec/CanonicalTests/Fire/Unit_Tests/LevelSetStartupChecks
                   ${INPUT_FILE} "${EXPECTED_MESSAGE}" "")
    set_tests_properties(${TEST_NAME} PROPERTIES LABELS "regression;fire")
endfunction(add_test_fire_levelset_startup)
add_test_fire_levelset_startup(FireLevelSetReinitClampBadValue_abort inputs_reinit_clamp_bad_value
    "erf.fire.levelset.reinit_clamp = \"step_begin\" is not one of")
add_test_fire_levelset_startup(FireLevelSetEpsViscUnstable_abort inputs_eps_visc_unstable
    "The level-set substep is unstable: erf.fire.levelset.eps_visc = 20 m")
add_test_fire_levelset_startup(FireLevelSetEpsViscFrontUnstable_abort inputs_eps_visc_front_unstable
    "The level-set substep is unstable: erf.fire.levelset.eps_visc_front = 20 m")
add_test_fire_levelset_startup(FireLevelSetCflUnstable_abort inputs_cfl_unstable_without_viscosity
    "No viscosity is stable at this cfl")
add_test_fire_levelset_startup(FireLevelSetEpsViscDiagonal_abort inputs_eps_visc_diagonal_unstable
    "The level-set substep is unstable: erf.fire.levelset.eps_visc = 12 m")
add_test_fire_levelset_startup(FireLevelSetEpsViscFrontBlend_abort inputs_eps_visc_front_blend_unstable
    "The level-set substep is unstable: the viscosity 13.53 m, blended from erf.fire.levelset.eps_visc_front = 14.75 m")
add_test_fire_levelset_startup(FireLevelSetEpsViscFrontWall_abort inputs_eps_visc_front_wall_first_order
    "The level-set substep is unstable: erf.fire.levelset.eps_visc_front = 12 m near the front, inside the WENO band, next to a non-burnable cell")
# erf.fire.levelset.reinit_clamp, old against new on one level-set deck
# (Exec/RegTests/FireReinitClamp): the "post_advection" front runs ahead of the
# "step_start" one, and with reinit_iters = 0 the two clamps give the same field
if(ERF_ENABLE_MPI AND NOT WIN32)
    add_test_fire_script(FireReinitClamp FireReinitClamp run_reinit_clamp.sh NRANKS 1)
endif()
# the checker's own pass/fail logic on synthetic fields. Pure Python.
add_test(FireReinitClamp_checker_SelfTest ${ERF_RANS_PYTHON}
    ${PROJECT_SOURCE_DIR}/Exec/RegTests/FireReinitClamp/check_reinit_clamp.py --self-test)
set_tests_properties(FireReinitClamp_checker_SelfTest
    PROPERTIES
    TIMEOUT 60
    PROCESSORS 1
    LABELS "unit;fire")
# every documented fire/dust key is read, every read key is documented, no deck sets an unread key.
# Labelled unit as well as docs: no CI job runs "ctest -L docs", so on the docs label alone this
# test went red for a whole PR without anyone seeing it (erf.fire.custom_fuel.*). It is pure
# Python and needs no binary, like the other checker self-tests under the unit label.
add_test(FireDustInputsDocs ${ERF_RANS_PYTHON} ${PROJECT_SOURCE_DIR}/Tests/check_fire_dust_inputs.py ${PROJECT_SOURCE_DIR})
set_tests_properties(FireDustInputsDocs PROPERTIES LABELS "docs;fire;unit" TIMEOUT 120)
# That checker's own pass/fail logic: an unread key is caught by nothing else, so a
# checker that stopped reporting one would pass the build in silence. Decks with a
# known verdict, including the key families whose names the code builds at run time
# (erf.fire.custom_fuel.<code>.*, erf.fire.firebreak.<n>.*), which used to read as
# dead. Pure Python, no ERF run.
add_test(FireDustInputsDocs_SelfTest ${ERF_RANS_PYTHON}
    ${PROJECT_SOURCE_DIR}/Tests/test_check_fire_dust_inputs.py ${PROJECT_SOURCE_DIR})
set_tests_properties(FireDustInputsDocs_SelfTest
    PROPERTIES
    TIMEOUT 300
    PROCESSORS 1
    LABELS "docs;fire;unit")
# Exec/Tools/farsite_winds_convert.py (erf.fire.farsite_winds grids to VTK,
# shapefile, KMZ, PDF, GeoPDF): a synthetic grid with known winds round-trips
# through VTK always and through each other format whose library is installed.
add_test(FarsiteWindsConvert_SelfTest ${ERF_RANS_PYTHON}
    ${PROJECT_SOURCE_DIR}/Exec/Tools/farsite_winds_convert.py --self-test)
set_tests_properties(FarsiteWindsConvert_SelfTest PROPERTIES TIMEOUT 120 PROCESSORS 1 LABELS "fire;unit")
# the mass-budget checker of FireDustMassConservation on a synthetic closed
# budget (passes) and a 10 % leak (must fail the budget check). Pure Python.
add_test(FireDustMassConservation_SelfTest ${ERF_RANS_PYTHON}
    ${PROJECT_SOURCE_DIR}/Exec/CanonicalTests/Hazard/FireDustMassConservation/check_mass_conservation.py --self-test)
set_tests_properties(FireDustMassConservation_SelfTest
    PROPERTIES
    TIMEOUT 60
    PROCESSORS 1
    LABELS "fire;unit")
if(ERF_ENABLE_DUST)
add_test_fire(FireRestart_dust_straight     FireRestart           inputs_dust_straight       40 NRANKS 1)
# the three fire-dust couplings applied once per step, in the right order
add_test_fire_check(FireDustCoupling_check  FireDustCoupling      inputs                     40 check_firedust.py NRANKS 1)
# each step's dust_diag.dat row written once across the final step and two restarts
if(ERF_ENABLE_MPI AND NOT WIN32)
add_test_fire_script(FireRestart_dust_rows  FireRestart           run_dust_rows.sh NRANKS 1)
# one rank and one box against two ranks and four boxes: the same dust fields,
# dust_diag.dat, CM budget and receptor sample (every dust CTest ran on one box
# until October 2026, so the ParallelCopy, average_down, box-edge and reduction
# paths had no coverage)
add_test_fire_script(FireDustCoupling_parity FireDustCoupling     run_dust_parity.sh NRANKS 2)
endif()
# dust inputs the kernels cannot use stop at start-up
# a dust start-up warning is printed once: DustLayer's copy of the dust inputs used
# to read them again, so every DustParams warning appeared twice
add_test_fire_abort(DustWarningOnce           FireRestart           inputs_dust_straight
    "erf.dust.threshold_A_coeff is read only with erf.dust.threshold_model = bagnold" "erf.dust.threshold_A_coeff=0.1")
set_tests_properties(DustWarningOnce PROPERTIES
    FAIL_REGULAR_EXPRESSION "threshold_A_coeff is read only.*threshold_A_coeff is read only")
add_test_fire_abort(DustBadBins_abort         FireRestart           inputs_dust_straight
    "erf.dust.n_size_bins must be >= 1" "erf.dust.n_size_bins=0")
# the same chop on the dust grid, with the fire off so the dust check is the one reached
add_test_fire_abort(DustBoxRatio_abort        FireRestart           inputs_dust_straight
    "DUST. Box sizes not divisible by grid_ratio: box 0 of level 0 is 7 x 7 cells.*Set amr.blocking_factor_x"
    "erf.fire.enable=false amr.max_grid_size_x=8 amr.max_grid_size_y=8")
add_test_fire_abort(DustZrefMismatch_abort    FireRestart           inputs_dust_straight
    "must equal erf.most.zref" "erf.most.zref=12.0")
# October 2026 audit: a missing road file used to hang every rank but the IO rank in
# a broadcast; the state carries one dust scalar, so separately transported bins
# overwrote the smoke and vapour slots; zref <= z0 silently gave u* = 0.4 |U|
add_test_fire_abort(DustRoadFileMissing_abort FireRestart           inputs_dust_straight
    "erf.dust.road_schedule_file cannot be opened" "erf.dust.road_schedule_file=missing_roads.csv")
add_test_fire_abort(DustTransportBins_abort   FireRestart           inputs_dust_straight
    "transport_bins_separately = true needs one state scalar per bin" "erf.dust.transport_bins_separately=true")
add_test_fire_abort(FireDustWindZref_abort    FireRestart           inputs_dust_straight
    "erf.fire_dust_wind_zref .* must exceed" "erf.fire_dust_wind_z0=10.0")   # zref follows erf.fire.wind_ref_ht = 6.1
# October 2026 validation: every new selector and range check, and the readers
# that used to warn and run on with the uniform value, proved to stop
add_test_fire_abort(DustThresholdModel_abort  FireRestart           inputs_dust_straight
    "erf.dust.threshold_model must be shao_lu or bagnold" "erf.dust.threshold_model=owen")
add_test_fire_abort(DustSaltationDiam_abort   FireRestart           inputs_dust_straight
    "erf.dust.saltation_diameter must be > 0" "erf.dust.saltation_diameter=0.0")
add_test_fire_abort(DustShaoLuAN_abort        FireRestart           inputs_dust_straight
    "erf.dust.shao_lu_A_N must be > 0" "erf.dust.shao_lu_A_N=0.0")
add_test_fire_abort(DustShaoLuGamma_abort     FireRestart           inputs_dust_straight
    "erf.dust.shao_lu_gamma must be >= 0" "erf.dust.shao_lu_gamma=-1.0e-4")
add_test_fire_abort(DustTerrainUstar_abort    FireRestart           inputs_dust_straight
    "erf.dust.terrain_ustar must be scale or loglaw" "erf.dust.terrain_ustar=log")
add_test_fire_abort(DustLumpedSettling_abort  FireRestart           inputs_dust_straight
    "erf.dust.lumped_settling must be mean or bin0" "erf.dust.lumped_settling=max")
add_test_fire_abort(DustAveraging_abort       FireRestart           inputs_dust_straight
    "erf.dust.averaging must be window or exponential" "erf.dust.averaging=running")
add_test_fire_abort(DustTooManyBins_abort     FireRestart           inputs_dust_straight
    "erf.dust.n_size_bins must be <= 8" "erf.dust.n_size_bins=9")
add_test_fire_abort(DustCmBudgetInt_abort     FireRestart           inputs_dust_straight
    "erf.dust.cm_budget_int must be >= 1" "erf.dust.cm_budget_int=0")
add_test_fire_abort(DustVisibilityKext_abort  FireRestart           inputs_dust_straight
    "erf.dust.visibility_k_ext must be > 0" "erf.dust.visibility_k_ext=0.0")
add_test_fire_abort(DustStelAveraging_abort   FireRestart           inputs_dust_straight
    "erf.dust.stel_averaging_s must be > 0" "erf.dust.stel_averaging_s=0.0")
add_test_fire_abort(DustMetalVarRemoved_abort FireRestart           inputs_dust_straight
    "erf.dust.phreeqc_metal_var was removed" "erf.dust.phreeqc_metal_var=metal_fraction")
add_test_fire_abort(DustMrfScRemoved_abort FireRestart              inputs_dust_straight
    "erf.dust_mrf_Sc_t was removed" "erf.dust_mrf_Sc_t=2.0")
add_test_fire_abort(DustTerrainFileMissing_abort FireRestart        inputs_dust_straight
    "erf.dust.terrain_file cannot be opened" "erf.dust.terrain_file=missing_terrain.asc")
add_test_fire_abort(DustRasterMissing_abort   FireRestart           inputs_dust_straight
    "surface raster cannot be opened" "erf.dust.crust_index_file=missing_crust.asc")
add_test_fire_abort(DustPhreeqcMissing_abort  FireRestart           inputs_dust_straight
    "PHREEQC file cannot be opened" "erf.dust.phreeqc_output_file=missing_phreeqc.csv")
# the haul-road mass is the AP-42 PM-10 factor on bin 0; the fire-wind height follows erf.fire.wind_ref_ht
add_test_fire_abort(DustRoadBinZero_abort      FireRestart           inputs_dust_straight
    "road mass into bin 0, whose diameter" "erf.dust.bin_diameters=50.0e-6\ 2.5e-6\ 7.0e-6")
add_test_fire_abort(FireDustWindZrefMatch_abort FireRestart          inputs_dust_straight
    "must match 6.1.* m, the height of the wind the fire hands over" "erf.fire_dust_wind_zref=10.0")
# with use_per_fuel_wind_ht the fire hands over the 6.096 m wind (until 2026-10
# the dust read it as the wind at wind_ref_ht, so zref = 6.1 ran here)
add_test_fire_abort(FireDustWindZrefPerFuel_abort FireRestart        inputs_dust_straight
    "must match 6.096.* m, the height of the wind the fire hands over" "erf.fire.use_per_fuel_wind_ht=true erf.fire_dust_wind_zref=6.1")
endif()
endif()

#=============================================================================
# RANS time-step limits (Exec/CanonicalTests/Canonical_RANS/Timestep_Limits)
#=============================================================================
# Largest stable time step of one closure under explicit anelastic, implicit
# anelastic and implicit compressible integration (Timestep_Limits). The
# sweep spins the closure up and climbs a ladder of steps, about 30 short ERF
# runs and 60-90 s on one rank in Release, so it is labelled dt_sweep and kept
# out of the regression label that the CI runs in Debug.
function(add_test_rans_dt TEST_NAME CLOSURE)
    set(_rans_root ${PROJECT_SOURCE_DIR}/Exec/CanonicalTests/Canonical_RANS)
    set(CURRENT_TEST_SOURCE_DIR ${_rans_root}/Timestep_Limits)
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_NAME})
    file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR})
    file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*")
    file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/")
    file(GLOB _rans_py "${ERF_SOURCE_DIR}/Exec/CanonicalTests/Canonical_RANS/*.py")
    file(COPY ${_rans_py} DESTINATION "${CURRENT_TEST_BINARY_DIR}/")
    file(COPY ${ERF_SOURCE_DIR}/Exec/Tools/erf_plotfile.py DESTINATION "${CURRENT_TEST_BINARY_DIR}/")

    # a 4 x 4 x 200 column: one rank
    if(ERF_ENABLE_MPI)
        # unquoted on purpose: the command is split by shlex inside \"...\", which nested
        # quotes would break; the test is registered for NOT WIN32 only
        set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}")
    else()
        unset(MPI_COMMANDS)
    endif()

    resolve_test_exe("" "${ERF_TEST_EXE}" TEST_EXE)

    set(test_log "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log")
    # The sweep's exit code is the verdict; its tables are echoed into the
    # ctest output so a failure shows every rung.
    set(test_command sh -c "rm -f ${CURRENT_TEST_BINARY_DIR}/SWEEP_FAILED && ( ${ERF_RANS_PYTHON} ${CURRENT_TEST_BINARY_DIR}/sweep_dt.py --exe ${TEST_EXE} --closure ${CLOSURE} --mpi-cmd \"${MPI_COMMANDS}\" --workdir ${CURRENT_TEST_BINARY_DIR}/dt_runs > ${test_log} 2>&1 || touch ${CURRENT_TEST_BINARY_DIR}/SWEEP_FAILED ) && cat ${test_log} && test ! -f ${CURRENT_TEST_BINARY_DIR}/SWEEP_FAILED")

    add_test(${TEST_NAME} ${test_command})
    set_tests_properties(${TEST_NAME}
        PROPERTIES
        TIMEOUT 3600
        PROCESSORS 1
        WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
        LABELS "rans;dt_sweep"
        ATTACHED_FILES_ON_FAIL "${test_log}"
    )
endfunction(add_test_rans_dt)

# The sweep's spin-up deck uses the FFT projection (the limits in its README
# were measured with it), so the entries exist only in FFT builds; the command
# runs through sh, so not on Windows.
if(ERF_ENABLE_FFT AND NOT WIN32)
    add_test_rans_dt(RANS_Timestep_Limits_kEqn      kEqn)
    add_test_rans_dt(RANS_Timestep_Limits_Deardorff Deardorff)
    add_test_rans_dt(RANS_Timestep_Limits_MRF       MRF)
endif()

# The fire coupling (ERF_ENABLE_FIRE): an unknown erf.fire.engine stops the run, the
# coupling inputs stop the run naming the key, and a build without the coupling
# refuses a fire deck instead of running it as plain ERF. The checks run in the
# constructor, so any deck serves (ERF's Terrain_Stress_Ridge).
if(ERF_ENABLE_MPI AND NOT WIN32)
set(_fire_abort_dir ${ERF_SOURCE_DIR}/Tests/test_files/Terrain_Stress_Ridge)
if(ERF_ENABLE_FIRE)
# flare always compiles its engine (FireEngine/) in, so the "no engine built in"
# abort cannot happen here; an unknown name can.
add_test_abort(Fire_Abort_EngineNotCompiledIn ${_fire_abort_dir} Terrain_Stress_Ridge.i
               "erf.fire.engine = .no_such_engine. is not one of"
               "erf.fire.enable=true erf.fire.engine=no_such_engine")
add_test_abort(Fire_Abort_OldPlotKey ${_fire_abort_dir} Terrain_Stress_Ridge.i
               "erf.fire_plot_int is now erf.fire.plot_int"
               "erf.fire_plot_int=10")
add_test_abort(Fire_Abort_CouplingType ${_fire_abort_dir} Terrain_Stress_Ridge.i
               "erf.fire.coupling_type = .tight. is not one of .passive., .lagged., .synchronous."
               "erf.fire.coupling_type=tight")
add_test_abort(Fire_Abort_AnchorLevel ${_fire_abort_dir} Terrain_Stress_Ridge.i
               "erf.fire.anchor_level must be -1 .the finest level. or a level >= 0, got -2"
               "erf.fire.anchor_level=-2")
add_test_abort(Fire_Abort_HeatFluxAlfg ${_fire_abort_dir} Terrain_Stress_Ridge.i
               "erf.fire.heat_flux_alfg must be > 0"
               "erf.fire.heat_flux_alfg=0")
add_test_abort(Fire_Abort_SourceMode ${_fire_abort_dir} Terrain_Stress_Ridge.i
               "erf.fire.source_mode = .replace. is not one of .overwrite., .add."
               "erf.fire.source_mode=replace")
add_test_abort(Fire_Abort_GridRatio ${_fire_abort_dir} Terrain_Stress_Ridge.i
               "erf.fire.grid_ratio must be >= 1, got 0"
               "erf.fire.grid_ratio=0")
add_test_abort(Fire_Abort_SmokeKeyWithoutSmoke ${_fire_abort_dir} Terrain_Stress_Ridge.i
               "erf.fire.smoke_emission_factor needs erf.fire.smoke_enable = true"
               "erf.fire.smoke_emission_factor=0.05")
# ... but with the fire switched off explicitly the key only warns, and the run goes on
add_test_abort(Fire_Warn_SmokeKeyFireOff ${_fire_abort_dir} Terrain_Stress_Ridge.i
               "erf.fire.smoke_emission_factor is ignored with erf.fire.enable = false"
               "erf.fire.enable=false erf.fire.smoke_emission_factor=0.05")
# A checkpoint that records no scalar layout predates the smoke scalar (one
# passive scalar); a build with the fire coupling (two) must refuse it rather
# than read every component after the scalars one slot off. The test writes a
# checkpoint, removes its ScalarLayout file and restarts from it.
set(_fl_dir ${CMAKE_CURRENT_BINARY_DIR}/test_files/Fire_Abort_CheckpointWithoutLayout)
file(MAKE_DIRECTORY ${_fl_dir})
file(GLOB _fl_files "${_fire_abort_dir}/*")
file(COPY ${_fl_files} DESTINATION ${_fl_dir})
resolve_test_exe("" "${ERF_TEST_EXE}" _fl_exe)
set(_fl_mpi "${ERF_SH_MPI_LAUNCHER}${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}")
add_test(NAME Fire_Abort_CheckpointWithoutLayout
         COMMAND sh -c "rm -rf chk_fl* && ${_fl_mpi} ${_fl_exe} Terrain_Stress_Ridge.i max_step=1 erf.plot_int_1=-1 erf.check_int=1 erf.check_file=chk_fl amrex.call_addr2line=0 > write.log 2>&1 && test -f chk_fl00001/ScalarLayout && rm chk_fl00001/ScalarLayout && ${_fl_mpi} ${_fl_exe} Terrain_Stress_Ridge.i max_step=2 erf.plot_int_1=-1 erf.check_int=-1 erf.restart=chk_fl00001 amrex.call_addr2line=0 2>&1 | tee Fire_Abort_CheckpointWithoutLayout.log"
         WORKING_DIRECTORY ${_fl_dir})
set_tests_properties(Fire_Abort_CheckpointWithoutLayout PROPERTIES
    LABELS "regression" TIMEOUT 600 PROCESSORS 1
    PASS_REGULAR_EXPRESSION "it records no scalar layout, so it predates the fire coupling's smoke scalar"
    ATTACHED_FILES_ON_FAIL "${_fl_dir}/write.log;${_fl_dir}/Fire_Abort_CheckpointWithoutLayout.log")
else()
add_test_abort(Fire_Abort_BuiltWithoutFire ${_fire_abort_dir} Terrain_Stress_Ridge.i
               "erf.fire.enable = true but this executable was built without the fire coupling"
               "erf.fire.enable=true")
endif()
endif()

# FLARE's structures (Structures/): the conductor model on the structure coupling
if(ERF_ENABLE_STRUCTURES)
  include(${PROJECT_SOURCE_DIR}/Structures/Tests/ERFStructuresTests.cmake)
  # the fragility ensemble driver's own tests (Python only)
  include(${PROJECT_SOURCE_DIR}/Structures/Tools/fragility_ensemble/FragilityEnsembleTests.cmake)
  # the fragility fit's own tests (Python only)
  include(${PROJECT_SOURCE_DIR}/Structures/Tools/fragility_fit/FragilityFitTests.cmake)
  # the broken-wire state at rest: the independent reference's own tests (Python only)
  include(${PROJECT_SOURCE_DIR}/Structures/Tools/broken_wire/BrokenWireTests.cmake)
endif()

# The dust coupling (ERF_ENABLE_DUST): flare always compiles its dust layer (Dust/), so an
# unknown erf.dust.model stops the run, and a build without the coupling refuses a dust
# deck instead of running it as plain ERF. The checks run at the end of start-up, so any
# deck serves (ERF's Terrain_Stress_Ridge).
if(ERF_ENABLE_MPI AND NOT WIN32)
set(_dust_abort_dir ${ERF_SOURCE_DIR}/Tests/test_files/Terrain_Stress_Ridge)
if(ERF_ENABLE_DUST)
add_test_abort(Dust_Abort_ModelNotCompiledIn ${_dust_abort_dir} Terrain_Stress_Ridge.i
               "erf.dust.model = .no_such_model. is not one of the dust models compiled in"
               "erf.dust.enable=true erf.dust.model=no_such_model")
else()
add_test_abort(Dust_Abort_BuiltWithoutDust ${_dust_abort_dir} Terrain_Stress_Ridge.i
               "erf.dust.enable = true but this executable was built without the dust coupling"
               "erf.dust.enable=true")
endif()
endif()

#=============================================================================
# The example models (Tests/FireEngineExample, Tests/DustModelExample and ERF's
# Tests/StructureModelExample) on ERF's model interfaces, in their own executable,
# flare_examples (Tests/CMakeLists.txt), so that flare keeps one fire engine and one dust
# model and a deck that names none gets FLARE's.
#=============================================================================
if(TARGET flare_examples AND ERF_ENABLE_MPI AND NOT WIN32)
set(_fex_src ${CMAKE_CURRENT_SOURCE_DIR}/FireEngineExample)
set(_dex_src ${CMAKE_CURRENT_SOURCE_DIR}/DustModelExample)
set(_fex_exe "$<TARGET_FILE:flare_examples>")
set(_dex_exe "$<TARGET_FILE:flare_examples>")
set(_smx_exe "$<TARGET_FILE:flare_examples>")
set(_fex_mpi "${ERF_SH_MPI_LAUNCHER}${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}")
set(_smx_mpi "${_fex_mpi}")

# The example engine: a coupled run writes the engine's fields into the ERF plotfile, and
# a restart from its mid-run checkpoint (engine state and the coupler's lagged fluxes)
# reproduces the final plotfile exactly.
if(FLARE_EXAMPLE_FIRE AND ERF_ENABLE_FCOMPARE)
set(_fex_dir ${CMAKE_CURRENT_BINARY_DIR}/test_files/Fire_ExampleEngine_Restart)
file(MAKE_DIRECTORY ${_fex_dir})
file(COPY ${_fex_src}/inputs_fire_example ${_fex_src}/input_sounding DESTINATION ${_fex_dir})
set(_fex_mpi "${ERF_SH_MPI_LAUNCHER}${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}")
add_test(NAME Fire_ExampleEngine_Restart
         COMMAND sh -c "rm -rf plt0* chk0* rst && ${_fex_mpi} ${_fex_exe} inputs_fire_example amrex.call_addr2line=0 > straight.log 2>&1 && grep -q fire_phi plt00008/Header && mkdir rst && cd rst && cp ../input_sounding . && ${_fex_mpi} ${_fex_exe} ../inputs_fire_example erf.restart=../chk00004 erf.check_int=-1 amrex.call_addr2line=0 > ../restart.log 2>&1 && cd .. && ${FCOMPARE_EXE} -r 0 --abs_tol 0 plt00008 rst/plt00008"
         WORKING_DIRECTORY ${_fex_dir})
set_tests_properties(Fire_ExampleEngine_Restart PROPERTIES
    LABELS "regression;fire" TIMEOUT 600 PROCESSORS 2
    ATTACHED_FILES_ON_FAIL "${_fex_dir}/straight.log;${_fex_dir}/restart.log")
endif()

# The example dust model: a constant flux from the ground into a closed box, so the dust
# ERF sums (RHO DUST) must be the mass the model emitted at every step; and a restart from
# the mid-run checkpoint (the model's state and the dust scalar) reproduces the final
# plotfile and the model's count exactly.
find_package(Python3 COMPONENTS Interpreter QUIET)
if(FLARE_EXAMPLE_DUST AND Python3_FOUND)
set(_dex_mpi "${ERF_SH_MPI_LAUNCHER}${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}")
set(_dex_mdir ${CMAKE_CURRENT_BINARY_DIR}/test_files/Dust_ExampleModel_Mass)
file(MAKE_DIRECTORY ${_dex_mdir})
file(COPY ${_dex_src}/inputs_dust_example ${_dex_src}/input_sounding DESTINATION ${_dex_mdir})
# 1e-4 kg/m2/s over the 2000 m x 2000 m ground
add_test(NAME Dust_ExampleModel_Mass
         COMMAND sh -c "rm -rf plt0* chk0* dust_example.csv && ${_dex_mpi} ${_dex_exe} inputs_dust_example erf.check_int=-1 amrex.call_addr2line=0 > run.log 2>&1 && grep -q rhoadv_dust plt00008/Header && \"${Python3_EXECUTABLE}\" ${_dex_src}/check_dust_example.py run.log dust_example.csv 1.0e-4 4.0e6"
         WORKING_DIRECTORY ${_dex_mdir})
set_tests_properties(Dust_ExampleModel_Mass PROPERTIES
    LABELS "regression;dust" TIMEOUT 600 PROCESSORS 2
    ATTACHED_FILES_ON_FAIL "${_dex_mdir}/run.log")
if(ERF_ENABLE_FCOMPARE)
set(_dex_rdir ${CMAKE_CURRENT_BINARY_DIR}/test_files/Dust_ExampleModel_Restart)
file(MAKE_DIRECTORY ${_dex_rdir})
file(COPY ${_dex_src}/inputs_dust_example ${_dex_src}/input_sounding DESTINATION ${_dex_rdir})
add_test(NAME Dust_ExampleModel_Restart
         COMMAND sh -c "rm -rf plt0* chk0* rst dust_example.csv && ${_dex_mpi} ${_dex_exe} inputs_dust_example amrex.call_addr2line=0 > straight.log 2>&1 && mkdir rst && cd rst && cp ../input_sounding . && ${_dex_mpi} ${_dex_exe} ../inputs_dust_example erf.restart=../chk00004 erf.check_int=-1 amrex.call_addr2line=0 > ../restart.log 2>&1 && cd .. && ${FCOMPARE_EXE} -r 0 --abs_tol 0 plt00008 rst/plt00008 && test \"$(tail -n 1 dust_example.csv)\" = \"$(tail -n 1 rst/dust_example.csv)\""
         WORKING_DIRECTORY ${_dex_rdir})
set_tests_properties(Dust_ExampleModel_Restart PROPERTIES
    LABELS "regression;dust" TIMEOUT 600 PROCESSORS 2
    ATTACHED_FILES_ON_FAIL "${_dex_rdir}/straight.log;${_dex_rdir}/restart.log")
endif()
endif()

# The fire view (StructureStep::fire) of ERF's example structure model, with the example
# fire engine
# Each run is the fire deck of Tests/FireEngineExample with the example model added.
# The engine burns a disk of radius 150 + 50 t m that releases 5e4 W/m2; a quarter of
# it is latent heat, which heat_release still counts. Each step the model prints the
# level and time of the fire's fields, the heat released and the area burned.
# check_fire_view.awk compares them with the disk to 3 %: fire cells of 10 m (5 m on
# level 1) resolve it to 1.4 % or better, and fields one fire step late or early miss
# by 11 % or more. Three layouts:
# - both on level 0: the fields belong to the fire step that ended at the model's step
#   time; a restart from the mid-run checkpoint must give the restart's time and the
#   heat of the step before it, which the example engine restores;
# - the fire on a level-1 patch, the model on level 0: as above;
# - the fire on level 0, the model on a level-1 patch: the fire is already at the end
#   of the coarse step;
# - the fire on a level-1 patch with three substeps of 0.1/3 s, the model on level 0:
#   the fire's time must equal the model's step time exactly, although three
#   substeps add up to 0.3 s where the coarse level has 0.30000000000000004 s.
if(FLARE_EXAMPLE_FIRE AND FLARE_EXAMPLE_STRUCTURES)
  set(_smf_common "erf.structures.model=example structure.example.force=50 0 0 structure.example.epsilon=50 fire.example.r0=150 fire.example.ros=50 fire.example.heat_flux=5.0e4 fire.example.nsub=1 fire.example.latent_fraction=0.25 erf.fire.grid_ratio=10 amr.max_grid_size=10 10 60 max_step=4 erf.check_int=-1 erf.plot_int_1=-1 amrex.call_addr2line=0")
  set(_smf_two "amr.max_level=1 erf.dt_ref_ratio=2 erf.refinement_indicators=box1 erf.box1.max_level=1")
  set(_smf_awk "awk -v r0=150 -v ros=50 -v q0=5.0e4 -v tol=0.03")
  set(_smf_names Structures_ExampleModel_FireView Structures_ExampleModel_FireViewFinerFire Structures_ExampleModel_FireViewFinerModel Structures_ExampleModel_FireViewExactTime)
  set(_smf_args_0 "structure.example.position=300 300 300")
  set(_smf_args_1 "${_smf_two} erf.box1.in_box_lo_indices_crse=5 5 0 erf.box1.in_box_hi_indices_crse=14 14 29 erf.fire.anchor_level=1 structure.example.level=0 structure.example.position=200 200 300")
  set(_smf_args_2 "${_smf_two} erf.box1.in_box_lo_indices_crse=0 0 0 erf.box1.in_box_hi_indices_crse=9 9 29 erf.fire.anchor_level=0 structure.example.level=1 structure.example.position=300 300 300")
  set(_smf_args_3 "amr.max_level=1 erf.dt_ref_ratio=3 erf.fixed_dt=0.1 erf.refinement_indicators=box1 erf.box1.max_level=1 erf.box1.in_box_lo_indices_crse=5 5 0 erf.box1.in_box_hi_indices_crse=14 14 29 erf.fire.anchor_level=1 structure.example.level=0 structure.example.position=200 200 300")
  set(_smf_check_0 "-v lev=0 -v dtf=0.5 -v lead=0")
  set(_smf_check_1 "-v lev=1 -v dtf=0.25 -v lead=0")
  set(_smf_check_2 "-v lev=0 -v dtf=0.5 -v lead=0.5")
  set(_smf_check_3 "-v lev=1 -v dtf=0.0333333333333333 -v lead=0")
  set(_smf_nlines_0 4)
  set(_smf_nlines_1 4)
  set(_smf_nlines_2 8)
  set(_smf_nlines_3 4)
  # The restart leg of the first run: checkpoint at step 2 (t = 1 s), restart, two steps
  set(_smf_more_0 "erf.check_int=2")
  set(_smf_rst_0 " && rm -rf rst && mkdir rst && cd rst && cp ../input_sounding . && ${_smx_mpi} ${_smx_exe} ../inputs_fire_example ${_smf_common} ${_smf_args_0} erf.restart=../chk00002 > ../restart.log 2>&1 && cd .. && ${_smf_awk} ${_smf_check_0} -v nlines=2 -f ${CMAKE_CURRENT_SOURCE_DIR}/StructureModelExample/check_fire_view.awk restart.log && grep -q 'structure example fire: t = 1 fire_lev = 0 fire_t = 1 ' restart.log")
  set(_smf_more_1 "")
  set(_smf_more_2 "")
  set(_smf_rst_1 "")
  set(_smf_rst_2 "")
  set(_smf_more_3 "")
  set(_smf_rst_3 "")
  foreach(_i RANGE 3)
    list(GET _smf_names ${_i} _smf_name)
    set(_smf_dir ${CMAKE_CURRENT_BINARY_DIR}/test_files/${_smf_name})
    file(MAKE_DIRECTORY ${_smf_dir})
    file(COPY ${_fex_src}/inputs_fire_example ${_fex_src}/input_sounding DESTINATION ${_smf_dir})
    add_test(NAME ${_smf_name}
             COMMAND sh -c "rm -rf chk0* && ${_smx_mpi} ${_smx_exe} inputs_fire_example ${_smf_common} ${_smf_args_${_i}} ${_smf_more_${_i}} > run.log 2>&1 && ${_smf_awk} ${_smf_check_${_i}} -v nlines=${_smf_nlines_${_i}} -f ${CMAKE_CURRENT_SOURCE_DIR}/StructureModelExample/check_fire_view.awk run.log${_smf_rst_${_i}}"
             WORKING_DIRECTORY ${_smf_dir})
    set_tests_properties(${_smf_name} PROPERTIES
        LABELS "regression;structures;fire" TIMEOUT 600 PROCESSORS 2
        ATTACHED_FILES_ON_FAIL "${_smf_dir}/run.log;${_smf_dir}/restart.log")
  endforeach()
endif()
endif()
