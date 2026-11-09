# FLARE's structures (Structures/): the conductor lines, towers and transformers on ERF's
# structure coupling (erf.structures.model = conductors). Tests/CTestList.cmake includes this
# file when ERF_ENABLE_STRUCTURES is ON. The cases are in Tests/test_files/Conductors_* and their
# golds in Tests/ERFGoldFiles, beside FLARE's fire cases; the scripts in Structures/Tests.

# The normwise column groups the conductor tests compare their logs with: the near-cancelling
# along-line pull on a tower must agree within its vector, and a real change of a large component
# or of a small cross-arm displacement must not. Pure CMake, no ERF run, hence the "unit" label.
add_test(CompareDataLogs_Groups ${CMAKE_COMMAND}
    "-DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/test_files/CompareDataLogs_Groups"
    -P ${PROJECT_SOURCE_DIR}/Structures/Tests/CompareDataLogsGroupsSelfTest.cmake)
set_tests_properties(CompareDataLogs_Groups
    PROPERTIES
    TIMEOUT 60
    PROCESSORS 1
    LABELS "unit;conductors")

# Power-line conductors (MoorDyn-C lines): single spans, circuits, towers and transformers,
# against golds and across restarts, and start-up checks of their inputs
#=============================================================================
if(ERF_ENABLE_STRUCTURES AND ERF_ENABLE_FFT AND ERF_ENABLE_MPI AND NOT WIN32 AND ERF_PRECISION STREQUAL "DOUBLE")
  # Run a conductor deck: the log LOG of one span (mid-span position, sag, offset, swing angle,
  # tensions, the wind handed to MoorDyn at the middle node) must match its gold GOLD_LOG to
  # SIGDIGITS digits, and the plotfile the gold PLOT_GOLD_NAME: the flow's gold, shared by the stub
  # and real-library tests, when nothing goes back into the flow, the case's own when the lines'
  # drag does (drag_on_flow). The bundled stub relaxes to the quasi-static angle while the real
  # MoorDyn swings about it, so each library has its own gold logs. TOTALS names the totals file
  # whose source integral is checked; EXTRA_LOGS lists further logs (space separated), each
  # compared with the gold of its file name followed by GOLD_SUFFIX (".gold" for the stub,
  # ".moordyn.gold" for the real library). OPTIONS (space-separated key=value inputs) follow the deck
  # on the command line: a variant of the deck, such as its old forms of a model against the golds
  # written before the model changed. RunConductors.cmake compares the components of one
  # vector, a tower's loads say, normwise (Tests/ConductorLogGroups.cmake).
  function(add_test_conductors TEST_NAME TEST_FILES_DIR PLTFILE LOG GOLD_LOG SIGDIGITS PLOT_GOLD_NAME)
      cmake_parse_arguments(ATC "" "TOTALS;EXTRA_LOGS;GOLD_SUFFIX;OPTIONS" "" ${ARGN})
      setup_test()
      set(PLOT_GOLD ${ERF_TEST_GOLD_FILES_DIRECTORY}/${PLOT_GOLD_NAME})
      resolve_test_exe("" "erf_exec" TEST_EXE)
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
          "-DPLOT_GOLD=${PLOT_GOLD}"
          "-DRTOL=${ERF_TEST_FCOMPARE_RTOL}"
          "-DATOL=${ERF_TEST_FCOMPARE_ATOL}"
          "-DLOG=${LOG}"
          "-DGOLD=${CURRENT_TEST_SOURCE_DIR}/${GOLD_LOG}"
          "-DSIGDIGITS=${SIGDIGITS}"
          "-DTOTALS=${ATC_TOTALS}"
          "-DEXTRA_LOGS=${ATC_EXTRA_LOGS}"
          "-DGOLD_DIR=${CURRENT_TEST_SOURCE_DIR}"
          "-DGOLD_SUFFIX=${ATC_GOLD_SUFFIX}"
          "-DOPTIONS=${ATC_OPTIONS}"
          -P ${PROJECT_SOURCE_DIR}/Structures/Tests/RunConductors.cmake)
      set_tests_properties(${TEST_NAME}
          PROPERTIES
          TIMEOUT 600
          PROCESSORS ${NP}
          WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
          LABELS "regression;conductors"
          ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/simulation.log")
  endfunction()
  # A deck whose towers move, run with the coupling steps ERF chooses and with SUBSTEPS: the columns
  # COLUMNS of the log LOG must agree row by row to 3 significant digits (2 units of the last).
  function(setup_test_coupling_convergence TEST_NAME TEST_FILES_DIR LOG COLUMNS SUBSTEPS)
      setup_test()
      resolve_test_exe("" "erf_exec" TEST_EXE)
      add_test(NAME ${TEST_NAME} COMMAND ${CMAKE_COMMAND}
          "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
          "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
          "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
          "-DNRANKS=${NP}"
          "-DTEST_EXE=${TEST_EXE}"
          "-DCONFIG=$<CONFIG>"
          "-DINPUT=${CURRENT_TEST_BINARY_DIR}/${TEST_FILES_DIR}.i"
          "-DWORKING_DIRECTORY=${CURRENT_TEST_BINARY_DIR}"
          "-DLOG=${LOG}"
          "-DCOLUMNS=${COLUMNS}"
          "-DSUBSTEPS=${SUBSTEPS}"
          "-DSIGDIGITS=3"
          "-DULPS=2"
          -P ${PROJECT_SOURCE_DIR}/Structures/Tests/RunCouplingConvergence.cmake)
      set_tests_properties(${TEST_NAME}
          PROPERTIES
          TIMEOUT 600
          PROCESSORS ${NP}
          WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
          LABELS "regression;conductors"
          ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/simulation_default.log;${CURRENT_TEST_BINARY_DIR}/simulation_fine.log")
  endfunction()
  if(ERF_MOORDYN_USE_STUB)
    # a prescribed 15 m/s crosswind handed to the line at every node, against the stub
    add_test_conductors(Conductors_PrescribedWind Conductors_PrescribedWind "plt00010" "S1.dat"
                        "S1.dat.gold" 8 Conductors_PrescribedWind)
    # the wind sampled from ERF's sheared crosswind (v = 15 m/s at the ground, 45 m/s at 300 m) at
    # the line's current nodes, against the stub; and the span's ASCE 74 design check (asce74.csv)
    add_test_conductors(Conductors_FlowWind Conductors_FlowWind "plt00010" "S1.dat"
                        "S1.dat.gold" 8 Conductors_FlowWind EXTRA_LOGS "conductors/asce74.csv" GOLD_SUFFIX ".gold")
    # the lines' drag put back into the flow: the source must integrate to the force on the air
    add_test_conductors(Conductors_DragOnFlow Conductors_DragOnFlow "plt00010" "S1.dat"
                        "S1.dat.gold" 8 Conductors_DragOnFlow TOTALS "conductors/total_load.dat")
  else()
    # the same two cases against the real MoorDyn-C: the span swings about the blowout angle;
    # compared to 6 digits, since the line integration is not bit-reproducible across compilers
    add_test_conductors(Conductors_FlowWind_MoorDyn Conductors_FlowWind "plt00010" "S1.dat"
                        "S1.dat.moordyn.gold" 6 Conductors_FlowWind EXTRA_LOGS "conductors/asce74.csv" GOLD_SUFFIX ".moordyn.gold")
    add_test_conductors(Conductors_DragOnFlow_MoorDyn Conductors_DragOnFlow "plt00010" "S1.dat"
                        "S1.dat.moordyn.gold" 6 Conductors_DragOnFlow_MoorDyn TOTALS "conductors/total_load.dat")
  endif()
  # The same span on a refined level: a static level-1 box that stops below the domain's top, its
  # steps subcycled; the lines step on it, take its wind and put their drag into its flow. The
  # source must integrate to the force on the air, the span's log and the two-level plotfile match
  # their golds, and the restart parity (its checkpointing run going on past the checkpoint) carries
  # the lines and their source across the checkpoint.
  if(ERF_MOORDYN_USE_STUB)
    add_test_conductors(Conductors_TwoLevel Conductors_TwoLevel "plt00010" "S1.dat"
                        "S1.dat.gold" 8 Conductors_TwoLevel TOTALS "conductors/total_load.dat")
    set(_two_level_restart Conductors_TwoLevel_Restart)
  else()
    add_test_conductors(Conductors_TwoLevel_MoorDyn Conductors_TwoLevel "plt00010" "S1.dat"
                        "S1.dat.moordyn.gold" 6 Conductors_TwoLevel_MoorDyn TOTALS "conductors/total_load.dat")
    set(_two_level_restart Conductors_TwoLevel_Restart_MoorDyn)
  endif()
  add_test_restart_parity(${_two_level_restart} Conductors_TwoLevel 5 10
      DATALOG "S1.dat conductors/total_load.dat" DATALOG_SIGDIGITS 10 CHK_LEG_END 10)
  set_tests_properties(${_two_level_restart} PROPERTIES LABELS "regression;restart-parity;conductors")
  # The capacities and the first-failure log on two levels (stub only: the pins are its still-air
  # and blown-out tensions). The span steps on level 1 (0.125 s, subcycled under level 0's 0.25 s),
  # so it is checked at every level-1 step and stamped with that step's end, never only at level 0's
  # (Ann's #5009 point for a model on a refined level). With rated_breaking_strength = 13892 N the
  # span's largest tension in S1.dat first passes it at t = 1.875 s, between two level-0 steps: the
  # first row of failures.dat must carry that time (a check at level 0's steps would give 2.0 s), the
  # case is required to cross off level 0's steps, and capacities.csv's peak demand must be the largest
  # tension in S1.dat. The pins sit close to the tension (13892 N is crossed at 13892.16 N after
  # 13891.89 N; 13890.45 N below, at 13890.53 N after 13890.39 N): a change of 1e-5 in the flow
  # can move a crossing, and then these capacities are set again from the new S1.dat.
  if(ERF_MOORDYN_USE_STUB)
    set(_cap_dir ${CMAKE_CURRENT_BINARY_DIR}/test_files/Conductors_Capacities_TwoLevel)
    file(MAKE_DIRECTORY ${_cap_dir})
    file(GLOB _cap_files ${CMAKE_CURRENT_SOURCE_DIR}/test_files/Conductors_TwoLevel/*)
    file(COPY ${_cap_files} DESTINATION ${_cap_dir})
    resolve_test_exe("" "erf_exec" _cap_exe)
    set(_cap_mpi "${ERF_SH_MPI_LAUNCHER}${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}")
    add_test(NAME Conductors_Capacities_TwoLevel
             COMMAND sh -c "rm -rf plt0* conductors S1.dat && ${_cap_mpi} ${_cap_exe} Conductors_TwoLevel.i conductors.S1.rated_breaking_strength=13892.0 amrex.call_addr2line=0 > run.log 2>&1 && awk -v rbs=13892.0 'FNR == 1 { f++ } f == 1 && FNR > 1 { if (t1 == \"\" && $10 > rbs) t1 = $1; if ($10 > peak) peak = $10 } f == 2 && FNR == 2 { tf = $1; kind = $2; name = $3 } f == 2 && FNR > 1 { rows++ } f == 3 && FNR == 2 { split($0, c, \",\"); pd = c[6] } END { coarse = (t1 / 0.25 == int(t1 / 0.25)); printf \"first crossing in S1.dat %s s, first row of failures.dat %s s (%s %s), %d row(s); peak tension %s N, capacities.csv peak demand %s N\\n\", t1, tf, kind, name, rows, peak, pd; exit !(t1 == 1.875 && tf == t1 && !coarse && rows == 1 && kind == \"conductor\" && name == \"S1\" && pd + 0 == peak + 0) }' S1.dat conductors/failures.dat conductors/capacities.csv"
             WORKING_DIRECTORY ${_cap_dir})
    set_tests_properties(Conductors_Capacities_TwoLevel PROPERTIES LABELS "regression;conductors;multilevel" TIMEOUT 600
                         PROCESSORS 2 ATTACHED_FILES_ON_FAIL "${_cap_dir}/run.log")
    # the records across a checkpoint at step 5 (t = 1.25 s): the conductor's first exceedance at 0.5 s
    # (rated_breaking_strength 13890.45 N) before it, the dead ends' at 1.875 s (hardware_rated_load
    # 13892 N) after it; the restarted run's failures.dat must be the straight run's
    add_test_restart_parity(Conductors_Capacities_Restart Conductors_TwoLevel 5 10
        COMMON_OPTIONS "conductors.S1.rated_breaking_strength=13890.45 conductors.S1.hardware_rated_load=13892.0"
        DATALOG "conductors/failures.dat conductors/capacities.csv S1.dat" DATALOG_SIGDIGITS 10 CHK_LEG_END 10)
    set_tests_properties(Conductors_Capacities_Restart PROPERTIES LABELS "regression;restart-parity;conductors;multilevel")
  endif()
  # A conductor span heated by a fire below it (conductors.thermal_model = ieee738; stub only, since the sag
  # check uses the stub's elastic parabola). FLARE's engine releases a prescribed 2e5 W/m^2 over a 40 m disc
  # under the middle of a 300 m span 30 m up, in still air at 300 K with passive coupling, so the conductor's
  # temperature follows its heat balance alone. check_fire_heating.py (README.md in the case folder)
  # recomputes from the deck's numbers, independently of the C++: the hottest node's irradiance at every step
  # from the fire cells and the node positions (1e-9), its temperature stepped with its own IEEE 738 terms
  # (1e-6 K over the logged 10 digits), the thermal strain of every row and the final sag from the stub's
  # parabola (1e-6). Three layouts of the fire and the lines on two levels, which also differ in which fire
  # step the lines see first (StructureStep::fire):
  # - one level: the first step sees no fire yet (the fire steps after the lines);
  # - the lines on a level-1 patch, the fire on level 0: the fire has stepped before the first fine step;
  # - the fire on a level-1 patch (5 m fire cells), the lines on level 0: the first step sees no fire yet.
  # The restart parity carries the temperatures, the spans' lengths and the fire across a checkpoint at 10 s.
  if(ERF_MOORDYN_USE_STUB AND ERF_ENABLE_FIRE)
    find_package(Python3 COMPONENTS Interpreter QUIET)
    if(Python3_Interpreter_FOUND)
      set(_fh_src ${CMAKE_CURRENT_SOURCE_DIR}/test_files/Conductors_FireHeating)
      set(_fh_two "amr.max_level=1 amr.ref_ratio_vect=2 2 2 erf.dt_ref_ratio=2 erf.coupling_type=TwoWay erf.refinement_indicators=box1 erf.box1.max_level=1 erf.box1.in_box_lo_indices_crse=10 6 0 erf.box1.in_box_hi_indices_crse=19 15 11")
      set(_fh_names Conductors_FireHeating Conductors_FireHeating_LinesOnFinerLevel Conductors_FireHeating_FireOnFinerLevel)
      set(_fh_opts_0 "")
      set(_fh_opts_1 "${_fh_two} erf.fire.anchor_level=0 conductors.anchor_level=-1")
      set(_fh_opts_2 "${_fh_two} erf.fire.anchor_level=1 conductors.anchor_level=0")
      set(_fh_check_0 "--dt 0.5 --fire-lags 1 --fire-dx 10")
      set(_fh_check_1 "--dt 0.25 --fire-lags 0 --fire-dx 10")
      set(_fh_check_2 "--dt 0.5 --fire-lags 1 --fire-dx 5")
      set(_fh_mpi "${ERF_SH_MPI_LAUNCHER}${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}")
      resolve_test_exe("" "erf_exec" _fh_exe)
      foreach(_i RANGE 2)
        list(GET _fh_names ${_i} _fh_name)
        set(_fh_dir ${CMAKE_CURRENT_BINARY_DIR}/test_files/${_fh_name})
        file(MAKE_DIRECTORY ${_fh_dir})
        file(COPY ${_fh_src}/Conductors_FireHeating.i ${_fh_src}/input_sounding DESTINATION ${_fh_dir})
        add_test(NAME ${_fh_name}
                 COMMAND sh -c "rm -rf conductors S1*.dat && ${_fh_mpi} ${_fh_exe} Conductors_FireHeating.i ${_fh_opts_${_i}} amrex.call_addr2line=0 > run.log 2>&1 && \"${Python3_EXECUTABLE}\" -I ${_fh_src}/check_fire_heating.py ${_fh_check_${_i}}"
                 WORKING_DIRECTORY ${_fh_dir})
        set_tests_properties(${_fh_name} PROPERTIES LABELS "regression;conductors;fire;multilevel" TIMEOUT 600
                             PROCESSORS 2 ATTACHED_FILES_ON_FAIL "${_fh_dir}/run.log")
      endforeach()
      add_test_restart_parity(Conductors_FireHeating_Restart Conductors_FireHeating 20 40
          COMMON_OPTIONS "erf.plot_file_1=plt erf.plot_vars_1=density theta x_velocity"
          DATALOG "conductors/conductor_temperature.dat S1.dat" DATALOG_SIGDIGITS 10 CHK_LEG_END 40)
      set_tests_properties(Conductors_FireHeating_Restart PROPERTIES LABELS "regression;restart-parity;conductors;fire")
    endif()
  endif()
  # The thermal model's start-up checks: an unknown model, a missing conductor property, a value out of range,
  # and a fire without the radiative fraction
  add_test_abort(Conductors_Thermal_UnknownModelAbort ${CMAKE_CURRENT_SOURCE_DIR}/test_files/Conductors_PrescribedWind
      Conductors_PrescribedWind.i "conductors.thermal_model must be none or ieee738, got 'cigre'" "conductors.thermal_model=cigre")
  add_test_abort(Conductors_Thermal_MissingPropertyAbort ${CMAKE_CURRENT_SOURCE_DIR}/test_files/Conductors_PrescribedWind
      Conductors_PrescribedWind.i "conductors.S1.heat_capacity is needed with conductors.thermal_model = ieee738"
      "conductors.thermal_model=ieee738 conductors.S1.absorptivity=0.5 conductors.S1.emissivity=0.5 conductors.S1.thermal_expansion=1.9e-5")
  add_test_abort(Conductors_Thermal_RangeAbort ${CMAKE_CURRENT_SOURCE_DIR}/test_files/Conductors_PrescribedWind
      Conductors_PrescribedWind.i "conductors.S1.absorptivity must be in .0, 1."
      "conductors.thermal_model=ieee738 conductors.S1.absorptivity=1.5 conductors.S1.emissivity=0.5 conductors.S1.heat_capacity=1310 conductors.S1.thermal_expansion=1.9e-5")
  if(ERF_ENABLE_FIRE AND ERF_MOORDYN_USE_STUB)
    # the fire deck without its radiative_fraction line (the command line cannot remove a key)
    set(_fn_dir ${CMAKE_CURRENT_BINARY_DIR}/test_files/Conductors_Thermal_FireNeedsFractionAbort)
    file(MAKE_DIRECTORY ${_fn_dir})
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_files/Conductors_FireHeating/input_sounding DESTINATION ${_fn_dir})
    resolve_test_exe("" "erf_exec" _fn_exe)
    add_test(NAME Conductors_Thermal_FireNeedsFractionAbort
             COMMAND sh -c "sed '/radiative_fraction/d' ${CMAKE_CURRENT_SOURCE_DIR}/test_files/Conductors_FireHeating/Conductors_FireHeating.i > nofraction.i && ${ERF_SH_MPI_LAUNCHER}${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ${_fn_exe} nofraction.i max_step=1 amrex.call_addr2line=0 > run.log 2>&1 && { echo 'the run did not stop'; exit 1; }; grep -q 'conductors.thermal.radiative_fraction is needed when the run has a fire' run.log"
             WORKING_DIRECTORY ${_fn_dir})
    set_tests_properties(Conductors_Thermal_FireNeedsFractionAbort PROPERTIES LABELS "regression;conductors;fire;abort" TIMEOUT 300
                         ATTACHED_FILES_ON_FAIL "${_fn_dir}/run.log")
  endif()
  # conductors.failure_stop stops the run at the first exceedance: the still-air span's ~13.5 kN
  # against a rated_breaking_strength of 10 kN, at the end of the first step (0.5 s)
  add_test_abort(Conductors_FailureStopAbort ${CMAKE_CURRENT_SOURCE_DIR}/test_files/Conductors_PrescribedWind
      Conductors_PrescribedWind.i "conductors.failure_stop: the first capacity exceedance, conductor S1 at t = 0.5"
      "conductors.S1.rated_breaking_strength=1.0e4 conductors.failure_stop=true")
  add_test_abort(Conductors_FailureStopNeedsCapacityAbort ${CMAKE_CURRENT_SOURCE_DIR}/test_files/Conductors_PrescribedWind
      Conductors_PrescribedWind.i "conductors.failure_stop = true needs a capacity to check" "conductors.failure_stop=true")
  set_tests_properties(Conductors_FailureStopAbort Conductors_FailureStopNeedsCapacityAbort PROPERTIES LABELS "regression;conductors")
  # conductors.ground_contact: the terrain is each line's MoorDyn seafloor. On flat ground, 330 m of
  # conductor over the 300 m chord sags 58 m from ends 30 m up and rests on the ground in the
  # crosswind; on the hills of Conductors_Terrain, both lines strung at 1.5 kN instead of 20 kN sag
  # onto the hillsides. Without contact both runs stop (a node below the ground leaves the domain;
  # Conductors_GroundContactOffAbort for the flat one). The restarted runs must carry the spans' ground records and their lines across the
  # checkpoint as the straight runs do; the lines reach the ground in the first steps
  if(ERF_MOORDYN_USE_STUB)
    set(_ground_flat Conductors_GroundContact_Restart)
    set(_ground_hills Conductors_GroundContact_Terrain_Restart)
  else()
    set(_ground_flat Conductors_GroundContact_Restart_MoorDyn)
    set(_ground_hills Conductors_GroundContact_Terrain_Restart_MoorDyn)
  endif()
  add_test_restart_parity(${_ground_flat} Conductors_PrescribedWind 5 10
      COMMON_OPTIONS "conductors.ground_contact=true conductors.S1.length=330"
      DATALOG "conductors/ground_contact.csv S1.dat" DATALOG_SIGDIGITS 10 CHK_LEG_END 10)
  add_test_restart_parity(${_ground_hills} Conductors_Terrain 5 10
      COMMON_OPTIONS "conductors.ground_contact=true conductors.L1.stringing_tension=1500 conductors.L2.stringing_tension=1500"
      DATALOG "conductors/ground_contact.csv conductors/L1_span2.dat" DATALOG_SIGDIGITS 10 CHK_LEG_END 10)
  set_tests_properties(${_ground_flat} ${_ground_hills} PROPERTIES LABELS "regression;restart-parity;conductors")
  # ground_contact with a MoorDyn step too long for the contact (five 60 m segments, moordyn_dt = 0.01 s)
  add_test_abort(Conductors_GroundContactStepAbort ${CMAKE_CURRENT_SOURCE_DIR}/test_files/Conductors_PrescribedWind
      Conductors_PrescribedWind.i "conductors.S1: with ground_contact MoorDyn's internal step"
      "conductors.ground_contact=true conductors.S1.segments=5 conductors.moordyn_dt=0.01")
  # the same flat-ground span without ground_contact: a node passes below the ground and leaves the domain
  add_test_abort(Conductors_GroundContactOffAbort ${CMAKE_CURRENT_SOURCE_DIR}/test_files/Conductors_PrescribedWind
      Conductors_PrescribedWind.i "has left the domain" "conductors.S1.length=330")
  set_tests_properties(Conductors_GroundContactStepAbort Conductors_GroundContactOffAbort PROPERTIES LABELS "regression;conductors")
  # A circuit: three phases 6 m apart hanging from insulator strings over two suspension towers,
  # and a shield wire clamped above them, across the same sheared crosswind. The middle phase's
  # middle span, its strings and their statistics, the shield wire's middle span and the
  # closest-approach statistics of the pairs P1-P2 and P2-SW and every span's ASCE 74 design check must
  # match their golds; the flow is Conductors_FlowWind's, since nothing goes back into it. Where along
  # two parallel lines they come closest is decided by millimetres, so the location is left to the
  # restart parity (one binary) and the unit tests, not to a gold.
  set(_circuit_logs "P2_insulators.dat P2_insulators_stats.csv SW_span2.dat conductors/separation_P1-P2_stats.csv conductors/separation_P2-SW_stats.csv")
  string(APPEND _circuit_logs " conductors/asce74.csv")
  set(_asce74_old_forms "conductors.asce74_wire_height=attachment conductors.asce74_inclined_spans=false")
  if(ERF_MOORDYN_USE_STUB)
    add_test_conductors(Conductors_Circuit Conductors_Circuit "plt00010" "P2_span2.dat"
                        "P2_span2.dat.gold" 8 Conductors_FlowWind EXTRA_LOGS "${_circuit_logs}" GOLD_SUFFIX ".gold")
    # the old forms of the ASCE 74 check (the attachment height, the whole weight on an inclined chord), with
    # the current elastic catenary, against their own asce74.csv golds
    add_test_conductors(Conductors_Circuit_Legacy Conductors_Circuit "plt00010" "conductors/asce74.csv"
                        "asce74.csv.legacy.gold" 8 Conductors_FlowWind OPTIONS "${_asce74_old_forms}")
    set(_circuit_restart Conductors_Circuit_Restart)
  else()
    add_test_conductors(Conductors_Circuit_MoorDyn Conductors_Circuit "plt00010" "P2_span2.dat"
                        "P2_span2.dat.moordyn.gold" 6 Conductors_FlowWind EXTRA_LOGS "${_circuit_logs}" GOLD_SUFFIX ".moordyn.gold")
    add_test_conductors(Conductors_Circuit_Legacy_MoorDyn Conductors_Circuit "plt00010" "conductors/asce74.csv"
                        "asce74.csv.legacy.moordyn.gold" 6 Conductors_FlowWind OPTIONS "${_asce74_old_forms}")
    set(_circuit_restart Conductors_Circuit_Restart_MoorDyn)
  endif()
  add_test_restart_parity(${_circuit_restart} Conductors_Circuit 5 10
      DATALOG "P2_span2.dat P2_insulators.dat SW_span2.dat conductors/separation.dat conductors/separation_P1-P2_stats.csv P2_insulators_stats.csv"
      DATALOG_SIGDIGITS 10)
  set_tests_properties(${_circuit_restart} PROPERTIES LABELS "regression;restart-parity;conductors")
  # Lines over hills, dead-ended on transformers: two sections from a hilltop transformer to two on
  # flat ground, on lattice towers, in a neutral boundary layer over terrain with the k-equation
  # RANS closure. The middle span of L1, the transformers' loads and clearances, the hilltop
  # transformer's statistics, the towers' drag and every span's gust factor and peak load from the
  # RANS k (gusts.csv) must match their golds; the flow is the run's own, since nothing goes back into it.
  set(_terrain_logs "conductors/transformers.dat conductors/transformer_T1_stats.csv conductors/towers.dat conductors/gusts.csv")
  if(ERF_MOORDYN_USE_STUB)
    add_test_conductors(Conductors_Terrain Conductors_Terrain "plt00010" "conductors/L1_span2.dat"
                        "L1_span2.dat.gold" 8 Conductors_Terrain EXTRA_LOGS "${_terrain_logs}" GOLD_SUFFIX ".gold")
    set(_terrain_restart Conductors_Terrain_Restart)
  else()
    # four digits: the lines start from MoorDyn's stationary initial-condition solve, whose
    # residual (the swing at time 0, 6.514e-4 deg) differs between compilers in the fifth digit,
    # and the RANS flow they then sample differs at roundoff
    add_test_conductors(Conductors_Terrain_MoorDyn Conductors_Terrain "plt00010" "conductors/L1_span2.dat"
                        "L1_span2.dat.moordyn.gold" 4 Conductors_Terrain EXTRA_LOGS "${_terrain_logs}" GOLD_SUFFIX ".moordyn.gold")
    set(_terrain_restart Conductors_Terrain_Restart_MoorDyn)
  endif()
  add_test_restart_parity(${_terrain_restart} Conductors_Terrain 5 10
      DATALOG "conductors/L1_span2.dat ${_terrain_logs}"
      DATALOG_SIGDIGITS 10)
  set_tests_properties(${_terrain_restart} PROPERTIES LABELS "regression;restart-parity;conductors")
  # The same lines in random gusts (gust_type = random) and in one travelling 1 - cos gust
  # (gust_type = event): every span and tower takes the RANS wind plus its gust. The middle span of
  # L1, the towers' loads and the gusts at each span's middle and each tower's top must match their
  # golds; the flow is Conductors_Terrain's, since nothing goes back into it. The restart parity
  # carries the random gusts' processes across the checkpoint; its checkpointing run goes on past the
  # checkpoint (CHK_LEG_END), so the restarted run must also drop the gust rows written after it.
  set(_gust_logs "conductors/towers.dat conductors/gust_series.dat")
  if(ERF_MOORDYN_USE_STUB)
    add_test_conductors(Conductors_RandomGusts Conductors_RandomGusts "plt00010" "conductors/L1_span2.dat"
                        "L1_span2.dat.gold" 8 Conductors_Terrain EXTRA_LOGS "${_gust_logs}" GOLD_SUFFIX ".gold")
    add_test_conductors(Conductors_EventGust Conductors_EventGust "plt00010" "conductors/L1_span2.dat"
                        "L1_span2.dat.gold" 8 Conductors_Terrain EXTRA_LOGS "${_gust_logs}" GOLD_SUFFIX ".gold")
    set(_gust_restart Conductors_RandomGusts_Restart)
  else()
    add_test_conductors(Conductors_RandomGusts_MoorDyn Conductors_RandomGusts "plt00010" "conductors/L1_span2.dat"
                        "L1_span2.dat.moordyn.gold" 4 Conductors_Terrain EXTRA_LOGS "${_gust_logs}" GOLD_SUFFIX ".moordyn.gold")
    add_test_conductors(Conductors_EventGust_MoorDyn Conductors_EventGust "plt00010" "conductors/L1_span2.dat"
                        "L1_span2.dat.moordyn.gold" 4 Conductors_Terrain EXTRA_LOGS "${_gust_logs}" GOLD_SUFFIX ".moordyn.gold")
    set(_gust_restart Conductors_RandomGusts_Restart_MoorDyn)
  endif()
  add_test_restart_parity(${_gust_restart} Conductors_RandomGusts 5 10
      DATALOG "conductors/L1_span2.dat ${_gust_logs}"
      DATALOG_SIGDIGITS 10 CHK_LEG_END 10)
  set_tests_properties(${_gust_restart} PROPERTIES LABELS "regression;restart-parity;conductors")
  # The same lines on towers that bend: each tower sways in its first mode under its members' drag
  # and its line's pull, and MoorDyn moves the cross-arms as coupled points. The towers' loads and
  # cross-arm displacements, the statistics of the hilltop line's first tower and the middle span of
  # L1 must match their golds; the flow is Conductors_Terrain's, since nothing goes back into it. The
  # restart parity carries the towers' sway and the coupling's log across the checkpoint; its
  # checkpointing run goes on past the checkpoint (CHK_LEG_END), so the restart must drop those rows.
  # The towers' old forms, for the legacy tests against the golds written before they changed.
  set(_tower_drag_old_forms "conductors.lattice.diagonal_wind_factor=false conductors.lattice.arm_outside_shaft=false")
  set(_moving_logs "conductors/towers.dat conductors/tower_L1_t1_stats.csv")
  if(ERF_MOORDYN_USE_STUB)
    add_test_conductors(Conductors_MovingTowers Conductors_MovingTowers "plt00010" "conductors/L1_span2.dat"
                        "L1_span2.dat.gold" 8 Conductors_Terrain EXTRA_LOGS "${_moving_logs}" GOLD_SUFFIX ".gold")
    # the old forms of the one-mode towers (no wind-direction factor, the cross-arm's drag through the shaft at
    # its height) give the towers' loads and sway and the span written before they changed
    add_test_conductors(Conductors_MovingTowers_Legacy Conductors_MovingTowers "plt00010" "conductors/towers.dat"
                        "towers.dat.legacy.gold" 8 Conductors_Terrain
                        EXTRA_LOGS "conductors/L1_span2.dat conductors/tower_L1_t1_stats.csv" GOLD_SUFFIX ".legacy.gold"
                        OPTIONS "${_tower_drag_old_forms}")
    set(_moving_restart Conductors_MovingTowers_Restart)
  else()
    add_test_conductors(Conductors_MovingTowers_MoorDyn Conductors_MovingTowers "plt00010" "conductors/L1_span2.dat"
                        "L1_span2.dat.moordyn.gold" 4 Conductors_Terrain EXTRA_LOGS "${_moving_logs}" GOLD_SUFFIX ".moordyn.gold")
    add_test_conductors(Conductors_MovingTowers_Legacy_MoorDyn Conductors_MovingTowers "plt00010" "conductors/towers.dat"
                        "towers.dat.legacy.moordyn.gold" 4 Conductors_Terrain
                        EXTRA_LOGS "conductors/L1_span2.dat conductors/tower_L1_t1_stats.csv" GOLD_SUFFIX ".legacy.moordyn.gold"
                        OPTIONS "${_tower_drag_old_forms}")
    set(_moving_restart Conductors_MovingTowers_Restart_MoorDyn)
  endif()
  # coupling.dat every 2 steps: the row after the step-5 checkpoint counts steps 5 and 6 across the restart
  add_test_restart_parity(${_moving_restart} Conductors_MovingTowers 5 10
      DATALOG "conductors/L1_span2.dat conductors/towers.dat conductors/tower_L1_t1_stats.csv conductors/coupling.dat"
      DATALOG_SIGDIGITS 10 CHK_LEG_END 10 COMMON_OPTIONS "conductors.diagnostics_int=2")
  set_tests_properties(${_moving_restart} PROPERTIES LABELS "regression;restart-parity;conductors")
  # The same lines on towers that bend as their frame model: every lattice tower stands on the
  # SubDyn frame of lattice_frame.dat (76 beam members on four fixed legs), its members' drag and
  # its line's pull carried onto the frame's nodes, and MoorDyn moves the cross-arms as coupled
  # points; lattice_members.dat gives the members' design data, so each is checked against its
  # strength. Case T has no joints where its diagonals cross, so its 9 m X diagonals are over the
  # bracing limit (KL/r 200), flagged, and govern on the elastic buckling branch at a utilisation
  # above 1: the gold tests the checks, not a sound design (Conductors_GeneratedTowers bolts its
  # diagonals where they cross). The towers' loads, footings (from the frame's support reactions), cross-arm
  # displacements and largest member utilisation, the statistics of L1's first tower, its members'
  # design strengths and the middle span of L1 must match their golds; the flow is
  # Conductors_Terrain's. The restart parity carries the frames' Newmark state, the members'
  # statistics and the coupling's log across the checkpoint, past which its checkpointing run goes on.
  set(_frame_logs "${_moving_logs} conductors/tower_L1_t1_members.csv")
  if(ERF_MOORDYN_USE_STUB)
    add_test_conductors(Conductors_FrameTowers Conductors_FrameTowers "plt00010" "conductors/L1_span2.dat"
                        "L1_span2.dat.gold" 8 Conductors_Terrain EXTRA_LOGS "${_frame_logs}" GOLD_SUFFIX ".gold")
    set(_frame_restart Conductors_FrameTowers_Restart)
  else()
    add_test_conductors(Conductors_FrameTowers_MoorDyn Conductors_FrameTowers "plt00010" "conductors/L1_span2.dat"
                        "L1_span2.dat.moordyn.gold" 4 Conductors_Terrain EXTRA_LOGS "${_frame_logs}" GOLD_SUFFIX ".moordyn.gold")
    set(_frame_restart Conductors_FrameTowers_Restart_MoorDyn)
  endif()
  if(NOT ERF_MOORDYN_USE_STUB)
    # The coupling of MoorDyn's lines and the frame towers must not depend on the coupling step: the
    # tensions of L1's middle span with the coupling steps ERF chooses (76 per ERF step here) and with
    # 160 agree to 3 digits. The insulator strings' axial modes (about 45 Hz) under the frames' stiff
    # cross-arms are not resolved by either count, so this holds only while the exchange between the
    # towers and the lines neither makes nor loses energy.
    setup_test_coupling_convergence(Conductors_FrameTowers_CouplingStep_MoorDyn Conductors_FrameTowers
                                    "conductors/L1_span2.dat" "tension_a tension_b max_tension" 160)
  endif()
  set(_frame_parity_logs "conductors/L1_span2.dat conductors/towers.dat conductors/tower_L1_t1_stats.csv")
  string(APPEND _frame_parity_logs " conductors/tower_L1_t1_members_stats.csv conductors/coupling.dat")
  add_test_restart_parity(${_frame_restart} Conductors_FrameTowers 5 10
      DATALOG "${_frame_parity_logs}" DATALOG_SIGDIGITS 10 CHK_LEG_END 10)
  set_tests_properties(${_frame_restart} PROPERTIES LABELS "regression;restart-parity;conductors")
  # The same lines on towers that bend as a frame ERF generates from the lattice type's dimensions
  # (8 panels of angles with crossed bracing and a cross-arm truss) with the steel at 500 C, every
  # member checked against its strength (ASCE 10-15) every step. The towers' loads, footings,
  # cross-arm displacements and largest member utilisation, the statistics of L1's first tower, its
  # members' design strengths and the middle span of L1 must match their golds; the flow is
  # Conductors_Terrain's. The restart parity carries the frames' Newmark state, the members'
  # statistics and the frames' log across the checkpoint; its checkpointing run goes on past the
  # checkpoint (CHK_LEG_END), so the restarted run must also drop the rows written after it.
  set(_generated_logs "conductors/towers.dat conductors/tower_L1_t1_stats.csv conductors/tower_L1_t1_members.csv")
  set(_tower_old_forms "${_tower_drag_old_forms} conductors.lattice.angle_principal_axes=false")
  if(ERF_MOORDYN_USE_STUB)
    add_test_conductors(Conductors_GeneratedTowers Conductors_GeneratedTowers "plt00010" "conductors/L1_span2.dat"
                        "L1_span2.dat.gold" 8 Conductors_Terrain EXTRA_LOGS "${_generated_logs}" GOLD_SUFFIX ".gold")
    # the old forms of the towers (no wind-direction factor, the cross-arm's drag through the shaft at its
    # height, every generated angle at MSpin 0) give the towers' loads and sway written before they changed
    add_test_conductors(Conductors_GeneratedTowers_Legacy Conductors_GeneratedTowers "plt00010" "conductors/towers.dat"
                        "towers.dat.legacy.gold" 8 Conductors_Terrain
                        EXTRA_LOGS "conductors/tower_L1_t1_stats.csv" GOLD_SUFFIX ".legacy.gold"
                        OPTIONS "${_tower_old_forms}")
    set(_generated_restart Conductors_GeneratedTowers_Restart)
  else()
    add_test_conductors(Conductors_GeneratedTowers_MoorDyn Conductors_GeneratedTowers "plt00010" "conductors/L1_span2.dat"
                        "L1_span2.dat.moordyn.gold" 4 Conductors_Terrain EXTRA_LOGS "${_generated_logs}" GOLD_SUFFIX ".moordyn.gold")
    add_test_conductors(Conductors_GeneratedTowers_Legacy_MoorDyn Conductors_GeneratedTowers "plt00010" "conductors/towers.dat"
                        "towers.dat.legacy.moordyn.gold" 4 Conductors_Terrain
                        EXTRA_LOGS "conductors/tower_L1_t1_stats.csv" GOLD_SUFFIX ".legacy.moordyn.gold"
                        OPTIONS "${_tower_old_forms}")
    set(_generated_restart Conductors_GeneratedTowers_Restart_MoorDyn)
  endif()
  set(_generated_parity_logs "conductors/L1_span2.dat conductors/towers.dat conductors/tower_L1_t1_stats.csv")
  string(APPEND _generated_parity_logs " conductors/tower_L1_t1_members_stats.csv conductors/tower_L1_t1_frame.dat")
  add_test_restart_parity(${_generated_restart} Conductors_GeneratedTowers 5 10
      DATALOG "${_generated_parity_logs}" DATALOG_SIGDIGITS 10 CHK_LEG_END 10)
  set_tests_properties(${_generated_restart} PROPERTIES LABELS "regression;restart-parity;conductors")
  # A circuit over hills that are an immersed boundary on a flat mesh: three phases and a shield wire
  # hang from two shared lattice towers that bend, standing on the hills' surface read from the
  # terrain file (the mesh's bottom is flat). The towers' loads and sway, the middle phase's middle
  # span and the transformers' loads must match their golds; the flow has its own gold, shared by
  # both libraries, since nothing goes back into it. The restart parity carries the shared towers'
  # sway and every line hanging from them across the checkpoint.
  set(_ib_logs "conductors/towers.dat conductors/transformers.dat")
  if(ERF_MOORDYN_USE_STUB)
    add_test_conductors(Conductors_ImmersedHills Conductors_ImmersedHills "plt00010" "conductors/L1b_span2.dat"
                        "L1b_span2.dat.gold" 8 Conductors_ImmersedHills EXTRA_LOGS "${_ib_logs}" GOLD_SUFFIX ".gold")
    set(_ib_restart Conductors_ImmersedHills_Restart)
  else()
    # four digits: the residual of MoorDyn's initial-condition solve and the LES flow differ between
    # compilers in the fifth digit; the towers' along-line pull at time 0, a few N left of span
    # pulls of 10 to 20 kN, is held to four digits of its tower's loads (ConductorLogGroups.cmake)
    add_test_conductors(Conductors_ImmersedHills_MoorDyn Conductors_ImmersedHills "plt00010" "conductors/L1b_span2.dat"
                        "L1b_span2.dat.moordyn.gold" 4 Conductors_ImmersedHills EXTRA_LOGS "${_ib_logs}" GOLD_SUFFIX ".moordyn.gold")
    set(_ib_restart Conductors_ImmersedHills_Restart_MoorDyn)
  endif()
  add_test_restart_parity(${_ib_restart} Conductors_ImmersedHills 5 10
      DATALOG "conductors/L1b_span2.dat conductors/L1a_span2.dat conductors/L1sw_span2.dat conductors/towers.dat conductors/transformers.dat"
      DATALOG_SIGDIGITS 10)
  set_tests_properties(${_ib_restart} PROPERTIES LABELS "regression;restart-parity;conductors")
  # Restart parity of the coupled case: ten steps straight, and five steps, a checkpoint and five
  # more from it. The plotfile (the flow and the lines' momentum source) and the line's logs, node
  # output, statistics and total load must come out the same; a line restarted from rest, or
  # statistics, step count or logs started over, all show. One test per library.
  if(ERF_MOORDYN_USE_STUB)
    set(_conductors_restart Conductors_Restart)
  else()
    set(_conductors_restart Conductors_Restart_MoorDyn)
  endif()
  add_test_restart_parity(${_conductors_restart} Conductors_DragOnFlow 5 10
      DATALOG "S1.dat S1_nodes.dat S1_stats.csv conductors/total_load.dat" DATALOG_SIGDIGITS 10)
  set_tests_properties(${_conductors_restart} PROPERTIES LABELS "regression;restart-parity;conductors")
  # The widths of the drag's kernel per line and per tower type: the coupled span gets a lattice tower at
  # its middle and a second span 200 m upwind, and each line and the tower type set the width the shared
  # conductors.epsilon gives in the other run (1.5 cells), while the shared one is 3 cells. The
  # plotfile (the flow and the lines' and tower's momentum source) must be the same to the bit; a line or
  # tower type whose own width were not read, or not used, would spread its drag over 3 cells.
  add_test_option_parity(Conductors_OwnWidthsParity Conductors_DragOnFlow "plt00010"
      COMMON_OPTIONS "ylo.dirichlet_file=../inflow_profile erf.input_sounding_file=../input_sounding conductors.lines=S1 S2 conductors.S1.towers=750. 500. 30. conductors.S1.length=151. 151. conductors.S1.tower_type=lattice conductors.tower_types=lattice conductors.lattice.base_width=6. conductors.lattice.top_width=1.5 conductors.lattice.solidity=0.2 conductors.lattice.arm_length=12. conductors.S2.end_a=600. 300. 30. conductors.S2.end_b=900. 300. 30. conductors.S2.length=301.5 conductors.S2.diameter=0.0281 conductors.S2.mass_per_length=1.628 conductors.S2.axial_stiffness=3.0e7 conductors.S2.output_root=S2"
      OFF_OPTIONS "conductors.epsilon=1.5"
      ON_OPTIONS  "conductors.epsilon=3.0 conductors.S1.epsilon=1.5 conductors.S2.epsilon=1.5 conductors.lattice.epsilon=1.5")
  set_tests_properties(Conductors_OwnWidthsParity PROPERTIES LABELS "regression;option-parity;conductors")
  # Start-up checks that stop a run before it starts, on the flow-sampled single span: the
  # conductor lines are listed in conductors.lines, and giving conductors.spans stops
  # the run with a message naming conductors.lines
  add_test_abort(Conductors_SpansKeyAbort
                 ${CMAKE_CURRENT_SOURCE_DIR}/test_files/Conductors_FlowWind
                 Conductors_FlowWind.i
                 "conductors.spans is not an input: list the conductor lines in conductors.lines"
                 "conductors.spans=S1")
  # an attachment outside the domain horizontally (x = 1600 m beyond the 1500 m domain) stops the
  # run when the line is placed on the ground, naming the line
  add_test_abort(Conductors_AttachmentOutsideAbort
                 ${CMAKE_CURRENT_SOURCE_DIR}/test_files/Conductors_FlowWind
                 Conductors_FlowWind.i
                 "conductors.S1.end_b at \\(1600[.0-9]*, 500[.0-9]*\\) lies outside the domain horizontally"
                 "conductors.S1.end_b=1600. 500. 30.")
  # a generated frame that cannot stand on its tower stops the run naming the tower and its cross-arm:
  # a cross-arm 40 m deep on a tower whose cross-arm is 30 m up
  add_test_abort(Conductors_GeneratedFrameAbort
                 ${CMAKE_CURRENT_SOURCE_DIR}/test_files/Conductors_GeneratedTowers
                 Conductors_GeneratedTowers.i
                 "tower L1_t1 \\(cross-arm 30[.0-9]* m\\): the cross-arm's depth \\(40 m\\) must be in"
                 "conductors.lattice.arm_depth=40")
  # a line's or a tower type's own drag width: refused without drag_on_flow, as conductors.epsilon is, and
  # refused when it is not positive, naming the key
  add_test_abort(Conductors_LineWidthNeedsDragAbort
                 ${CMAKE_CURRENT_SOURCE_DIR}/test_files/Conductors_FlowWind
                 Conductors_FlowWind.i
                 "conductors.S1.epsilon needs conductors.drag_on_flow = true"
                 "conductors.S1.epsilon=1.5")
  add_test_abort(Conductors_TowerWidthAbort
                 ${CMAKE_CURRENT_SOURCE_DIR}/test_files/Conductors_Terrain
                 Conductors_Terrain.i
                 "conductors.lattice.epsilon must be positive"
                 "conductors.drag_on_flow=true conductors.lattice.epsilon=0")
  # a member design file that cannot be read stops the run naming its key
  add_test_abort(Conductors_MemberFileAbort
                 ${CMAKE_CURRENT_SOURCE_DIR}/test_files/Conductors_FrameTowers
                 Conductors_FrameTowers.i
                 "conductors.lattice.member_file = missing.dat: cannot read the member design file 'missing.dat'"
                 "conductors.lattice.member_file=missing.dat")
  # gusts from the RANS k: on a run without the k-equation RANS (the laminar flow-sampled span), a gust
  # input given without a gust_type, gusts with a prescribed wind, random gusts with the drag put into
  # the flow, and an event without its time each stop the run, naming the input
  add_test_abort(Conductors_GustNeedsRANSAbort
                 ${CMAKE_CURRENT_SOURCE_DIR}/test_files/Conductors_FlowWind
                 Conductors_FlowWind.i
                 "conductors.gust_type = factor needs the k-equation RANS \\(erf.rans_type = kEqn\\) on the conductors' anchor level 0"
                 "conductors.gust_type=factor")
  add_test_abort(Conductors_GustKeyAbort
                 ${CMAKE_CURRENT_SOURCE_DIR}/test_files/Conductors_FlowWind
                 Conductors_FlowWind.i
                 "conductors.gust_peak_factor needs an conductors.gust_type"
                 "conductors.gust_peak_factor=3")
  add_test_abort(Conductors_GustPrescribedWindAbort
                 ${CMAKE_CURRENT_SOURCE_DIR}/test_files/Conductors_PrescribedWind
                 Conductors_PrescribedWind.i
                 "conductors.gust_type = factor takes k from the flow, so it cannot be used with conductors.prescribed_velocity"
                 "conductors.gust_type=factor")
  add_test_abort(Conductors_GustDragOnFlowAbort
                 ${CMAKE_CURRENT_SOURCE_DIR}/test_files/Conductors_Terrain
                 Conductors_Terrain.i
                 "conductors.gust_type = random adds its gust to the wind of the lines and towers only; it cannot be used"
                 "conductors.gust_type=random conductors.drag_on_flow=1")
  add_test_abort(Conductors_GustEventTimeAbort
                 ${CMAKE_CURRENT_SOURCE_DIR}/test_files/Conductors_Terrain
                 Conductors_Terrain.i
                 "conductors.gust_type = event needs conductors.gust_event_time"
                 "conductors.gust_type=event conductors.gust_event_speed=20")
  # amr.max_level = 1 with nothing tagged: the lines' level does not exist at the start
  add_test_abort(Conductors_AnchorLevelMissingAbort
                 ${CMAKE_CURRENT_SOURCE_DIR}/test_files/Conductors_TwoLevel
                 Conductors_TwoLevel.i
                 "which does not exist now .the finest level is 0."
                 "erf.box1.max_level=0")
  set_tests_properties(Conductors_SpansKeyAbort Conductors_AttachmentOutsideAbort Conductors_GeneratedFrameAbort
                       Conductors_MemberFileAbort Conductors_GustNeedsRANSAbort Conductors_GustKeyAbort
                       Conductors_GustPrescribedWindAbort Conductors_GustDragOnFlowAbort Conductors_GustEventTimeAbort
                       Conductors_AnchorLevelMissingAbort
                       PROPERTIES LABELS "regression;conductors")
endif()

#=============================================================================
# The conductor model on more than one level (erf-model/ERF#5009 review): the lines' drag must
# reach the composite flow whole. Conductors_MomentumBudget.i is a periodic box of uniform air
# with no gravity, Coriolis force or diffusion, blowing across one line at 10 m/s: nothing but
# the drag changes the air's y-momentum, so the change of the volume integral of y_velocity
# between plt00000 and plt00010 (density 1; fvolumesum integrates the composite, each place on
# its finest level) must equal the impulse the model reports, the sum over the lines' steps of
# source_y in conductors/total_load.dat times the step. With the drag on a level a finer level
# covers, averaging down removes it (Ann's #5009 point): with ERF's check switched off and level 1
# over the whole box, the lines on level 0 lost all of their drag. ERF's structure coupling now
# stops that run (Conductors::support), which Conductors_Abort_FinerLevel* check. ERF does not
# synchronise momentum fluxes at coarse/fine faces (no momentum reflux): the y-momentum crossing
# a face normal to y, mostly the pressure term, is not the same on the two levels, and with the
# finer patch's y-faces 75 to 150 m from the line the budget missed by 23 % to 111 %, a limit of
# the multi-level coupling, not of the model. The finer levels here have coarse/fine faces only
# normal to x or z, where the y-momentum flux has no pressure term.
#=============================================================================
if(ERF_ENABLE_STRUCTURES AND ERF_ENABLE_MPI AND ERF_ENABLE_FCOMPARE AND NOT WIN32 AND ERF_PRECISION STREQUAL "DOUBLE")
  get_filename_component(_cmb_tools "${FCOMPARE_EXE}" DIRECTORY)
  get_filename_component(_cmb_fcmp_name "${FCOMPARE_EXE}" NAME)
  string(REPLACE "fcompare" "fvolumesum" _cmb_fvs_name "${_cmb_fcmp_name}")
  set(_cmb_fvs "${_cmb_tools}/${_cmb_fvs_name}")
  resolve_test_exe("" "erf_exec" _cmb_exe)
  set(_cmb_mpi "${ERF_SH_MPI_LAUNCHER}${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}")
  set(_cmb_src ${CMAKE_CURRENT_SOURCE_DIR}/test_files/Conductors_MomentumBudget)
  set(_cmb_two "amr.max_level=1 amr.ref_ratio_vect=2 2 2 erf.dt_ref_ratio=2 erf.coupling_type=TwoWay erf.refinement_indicators=box1 erf.box1.max_level=1")
  # the line 200 m up in a box 300 m high, so that a finer level can lie below its drag's reach
  set(_cmb_tall "geometry.prob_hi=600 400 300 amr.n_cell=24 16 12")

  # The first step's source against the drag of the line at rest in the 10 m/s wind,
  # 1/2 rho Cd D L U^2 = 1/2 1 1 6 300 100 = 9e4 N, so that a wrong drag, which moves the applied
  # and the reported source together, fails too: exactly (to 1e-6) with the stub, whose line holds
  # its quasi-static shape; with MoorDyn-C the line starts to move downwind within the step, which
  # lowers the drag relative to it (6.18e4 N measured), so between half of it and all of it.
  if(ERF_MOORDYN_USE_STUB)
    set(_cmb_first_lo 89999.91)
    set(_cmb_first_hi 90000.09)
  else()
    set(_cmb_first_lo 45000.0)
    set(_cmb_first_hi 90000.09)
  endif()

  # run Conductors_MomentumBudget.i with OPTIONS and require the y-momentum change to match the
  # reported impulse to the relative tolerance RTOL_STUB (the stub) or RTOL_MOORDYN (MoorDyn-C).
  # The stub's line holds its quasi-static shape; MoorDyn-C's swings and bobs, and the velocities
  # it sends through the coarse/fine faces carry momentum fluxes the levels do not synchronise
  # (the same limit as above, its advective part), so its two-level budgets close less tightly.
  function(add_test_conductor_budget TEST_NAME RTOL_STUB RTOL_MOORDYN OPTIONS)
    if(ERF_MOORDYN_USE_STUB)
      set(RTOL ${RTOL_STUB})
    else()
      set(RTOL ${RTOL_MOORDYN})
    endif()
    set(_dir ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_NAME})
    file(MAKE_DIRECTORY ${_dir})
    file(COPY ${_cmb_src}/Conductors_MomentumBudget.i DESTINATION ${_dir})
    add_test(NAME ${TEST_NAME}
             COMMAND sh -c "rm -rf plt0* conductors S1* && ${_cmb_mpi} ${_cmb_exe} Conductors_MomentumBudget.i ${OPTIONS} amrex.call_addr2line=0 > run.log 2>&1 && a=$(${_cmb_fvs} -v y_velocity plt00000 | awk '/integral of y_velocity/ { print $NF }') && b=$(${_cmb_fvs} -v y_velocity plt00010 | awk '/integral of y_velocity/ { print $NF }') && test -n \"$a\" && test -n \"$b\" && awk -v a=\"$a\" -v b=\"$b\" -v tol=${RTOL} -v flo=${_cmb_first_lo} -v fhi=${_cmb_first_hi} 'NR == 2 { f = -$9 } NR > 1 { s += $9 * ($1 - tp); tp = $1; n++ } END { d = (b - a) - s; if (d < 0) d = -d; m = (s < 0) ? -s : s; printf \"first step drag %.7e N (allowed %g to %g; 1/2 rho Cd D L U^2 = 9e4), change of the y-momentum %.7e, impulse of the drag %.7e (%d steps), relative difference %.2e, tolerance %g\\n\", f, flo, fhi, b - a, s, n, d / m, tol; exit !(n >= 10 && f >= flo && f <= fhi && d <= tol * m) }' conductors/total_load.dat"
             WORKING_DIRECTORY ${_dir})
    set_tests_properties(${TEST_NAME} PROPERTIES LABELS "regression;conductors;multilevel" TIMEOUT 600 PROCESSORS 2
                         ATTACHED_FILES_ON_FAIL "${_dir}/run.log")
  endfunction()

  # one level: 5e-7 measured with the stub, 3.7e-7 with MoorDyn-C
  add_test_conductor_budget(Conductors_MomentumBudget 1.0e-5 1.0e-5 "")
  # the lines on the finest level (anchor_level -1), a slab over the whole box up to 150 m that
  # holds them: 1.3e-5 measured with the stub, the flux through the slab's top face; 1.0e-3 with
  # MoorDyn-C
  add_test_conductor_budget(Conductors_MomentumBudget_FinestLevel 1.0e-4 3.0e-3
      "${_cmb_two} ${_cmb_tall} erf.box1.in_box_lo_indices_crse=0 0 0 erf.box1.in_box_hi_indices_crse=23 15 5")
  # the lines on the finest level, a patch over the whole width in y between x = 75 and 525 m:
  # its side faces are normal to x, where the y-momentum flux has no pressure term; 1.9e-6 measured
  # with the stub, 5.3e-4 with MoorDyn-C
  add_test_conductor_budget(Conductors_MomentumBudget_FinestLevelPatch 1.0e-5 2.0e-3
      "${_cmb_two} erf.box1.in_box_lo_indices_crse=3 0 0 erf.box1.in_box_hi_indices_crse=20 15 5")
  # the lines on level 0, a finer slab below their drag's reach (0 to 50 m): 1.1e-6 measured with
  # the stub, 9.8e-5 with MoorDyn-C
  add_test_conductor_budget(Conductors_MomentumBudget_FinerLevelBelow 1.0e-5 5.0e-4
      "${_cmb_two} ${_cmb_tall} conductors.anchor_level=0 conductors.S1.end_a=150 200 200 conductors.S1.end_b=450 200 200 erf.box1.in_box_lo_indices_crse=0 0 0 erf.box1.in_box_hi_indices_crse=23 15 1")

  # the lines on level 0 and a finer level over them: the drag would be removed, so the run stops
  add_test_abort(Conductors_Abort_FinerLevelCoversLines ${_cmb_src} Conductors_MomentumBudget.i
      "structure model .conductors. couples to level 0, but level 1 covers part of its support"
      "${_cmb_two} conductors.anchor_level=0 erf.box1.in_box_lo_indices_crse=3 5 0 erf.box1.in_box_hi_indices_crse=20 10 5")
  # a finer level 100 m downwind of the line (y 300 to 325 m) that holds no node but the last
  # 12.5 m of the drag kernel's reach of 3 epsilon = 112.5 m: a support without the kernel's
  # reach, or with 2 epsilon of it, would miss it
  add_test_abort(Conductors_Abort_FinerLevelCoversDragTail ${_cmb_src} Conductors_MomentumBudget.i
      "structure model .conductors. couples to level 0, but level 1 covers part of its support"
      "${_cmb_two} conductors.anchor_level=0 erf.box1.in_box_lo_indices_crse=0 12 0 erf.box1.in_box_hi_indices_crse=23 12 5")
  # an embedded-boundary terrain: ERF hands the model no mesh heights or cell volumes there
  add_test_abort(Conductors_Abort_EBTerrain ${_cmb_src} Conductors_MomentumBudget.i
      "erf.structures.model = conductors does not support erf.terrain_type = EB"
      "erf.terrain_type=EB")
  set_tests_properties(Conductors_Abort_FinerLevelCoversLines Conductors_Abort_FinerLevelCoversDragTail
                       PROPERTIES LABELS "regression;conductors;multilevel")
  set_tests_properties(Conductors_Abort_EBTerrain PROPERTIES LABELS "regression;conductors")

  # lines without drag_on_flow on level 0 under a finer level: nothing goes into the flow, so the
  # run goes on and its flow is the flow without the model, to the last bit
  set(_cow_dir ${CMAKE_CURRENT_BINARY_DIR}/test_files/Conductors_OneWayUnderFinerLevel)
  file(MAKE_DIRECTORY ${_cow_dir})
  file(COPY ${_cmb_src}/Conductors_MomentumBudget.i DESTINATION ${_cow_dir})
  set(_cow_opts "${_cmb_two} conductors.anchor_level=0 conductors.drag_on_flow=false erf.box1.in_box_lo_indices_crse=3 5 0 erf.box1.in_box_hi_indices_crse=20 10 5 erf.plot_vars_1=density x_velocity y_velocity z_velocity amrex.call_addr2line=0")
  add_test(NAME Conductors_OneWayUnderFinerLevel
           COMMAND sh -c "rm -rf plt0* conductors S1* nomodel && ${_cmb_mpi} ${_cmb_exe} Conductors_MomentumBudget.i ${_cow_opts} > run.log 2>&1 && test -s S1.dat && mkdir nomodel && grep -v '^erf.structures.model' Conductors_MomentumBudget.i > nomodel/Conductors_MomentumBudget.i && cd nomodel && ${_cmb_mpi} ${_cmb_exe} Conductors_MomentumBudget.i ${_cow_opts} > run.log 2>&1 && cd .. && ${FCOMPARE_EXE} -r 0 --abs_tol 0 plt00010 nomodel/plt00010"
           WORKING_DIRECTORY ${_cow_dir})
  set_tests_properties(Conductors_OneWayUnderFinerLevel PROPERTIES LABELS "regression;conductors;multilevel" TIMEOUT 600 PROCESSORS 2
                       ATTACHED_FILES_ON_FAIL "${_cow_dir}/run.log;${_cow_dir}/nomodel/run.log")
endif()
