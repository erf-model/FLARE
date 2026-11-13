# The sources of FLARE's structures. FLARE's CMakeLists.txt (flare_add_models) includes
# this file with ERF_ENABLE_STRUCTURES and calls erf_structure_model_add with the target to
# compile them into; ERF_STRUCTURE_MODEL_DIR holds this directory during the call. The conductor model (Conductors/ERF_ConductorModel.cpp)
# registers on ERF's structure coupling as erf.structures.model = conductors.
function(erf_structure_model_add target)
  if(NOT ERF_ENABLE_STRUCTURES)
    message(FATAL_ERROR "FLARE's Structures/ directory is compiled in with ERF_ENABLE_STRUCTURES=ON "
                        "(it needs the MoorDyn set-up of Structures/ERFStructures.cmake)")
  endif()
  target_sources(${target} PRIVATE
    ${ERF_STRUCTURE_MODEL_DIR}/MoorDyn/ERF_MoorDynSystem.cpp
    ${ERF_STRUCTURE_MODEL_DIR}/Towers/ERF_MemberDrag.cpp
    ${ERF_STRUCTURE_MODEL_DIR}/Towers/ERF_TowerInputs.cpp
    ${ERF_STRUCTURE_MODEL_DIR}/Towers/ERF_Tower.cpp
    ${ERF_STRUCTURE_MODEL_DIR}/Towers/ERF_TowerDynamics.cpp
    ${ERF_STRUCTURE_MODEL_DIR}/Towers/ERF_FrameInputs.cpp
    ${ERF_STRUCTURE_MODEL_DIR}/Towers/ERF_Frame.cpp
    ${ERF_STRUCTURE_MODEL_DIR}/Towers/ERF_FrameDynamics.cpp
    ${ERF_STRUCTURE_MODEL_DIR}/Towers/ERF_FrameTower.cpp
    ${ERF_STRUCTURE_MODEL_DIR}/Towers/ERF_MemberChecks.cpp
    ${ERF_STRUCTURE_MODEL_DIR}/Towers/ERF_LatticeFrame.cpp
    ${ERF_STRUCTURE_MODEL_DIR}/Conductors/ERF_ConductorInputs.cpp
    ${ERF_STRUCTURE_MODEL_DIR}/Conductors/ERF_ASCE74.cpp
    ${ERF_STRUCTURE_MODEL_DIR}/Conductors/ERF_Gusts.cpp
    ${ERF_STRUCTURE_MODEL_DIR}/Conductors/ERF_MoorDynInputWriter.cpp
    ${ERF_STRUCTURE_MODEL_DIR}/Conductors/ERF_ConductorLine.cpp
    ${ERF_STRUCTURE_MODEL_DIR}/Conductors/ERF_Transformer.cpp
    ${ERF_STRUCTURE_MODEL_DIR}/Conductors/ERF_Capacities.cpp
    ${ERF_STRUCTURE_MODEL_DIR}/Conductors/ERF_Intensity.cpp
    ${ERF_STRUCTURE_MODEL_DIR}/Conductors/ERF_GroundContact.cpp
    ${ERF_STRUCTURE_MODEL_DIR}/Conductors/ERF_LineBreaks.cpp
    ${ERF_STRUCTURE_MODEL_DIR}/Conductors/ERF_TowerFailures.cpp
    ${ERF_STRUCTURE_MODEL_DIR}/Conductors/ERF_LineThermal.cpp
    ${ERF_STRUCTURE_MODEL_DIR}/Conductors/ERF_TowerThermal.cpp
    ${ERF_STRUCTURE_MODEL_DIR}/Conductors/ERF_Conductors.cpp
    ${ERF_STRUCTURE_MODEL_DIR}/Conductors/ERF_ConductorModel.cpp
  )
  target_include_directories(${target} PUBLIC
    $<BUILD_INTERFACE:${ERF_STRUCTURE_MODEL_DIR}/MoorDyn>
    $<BUILD_INTERFACE:${ERF_STRUCTURE_MODEL_DIR}/Towers>
    $<BUILD_INTERFACE:${ERF_STRUCTURE_MODEL_DIR}/Conductors>
    $<INSTALL_INTERFACE:include>)
  # read by ERF_MoorDynSystem.cpp only
  target_compile_definitions(${target} PRIVATE ERF_MOORDYN_VERSION="${ERF_MOORDYN_VERSION}")
  if(ERF_MOORDYN_USE_STUB)
    target_compile_definitions(${target} PRIVATE ERF_MOORDYN_USE_STUB)
    target_link_libraries(${target} PUBLIC erf_moordyn_stub)
  else()
    target_link_libraries(${target} PUBLIC MoorDyn::moordyn)
  endif()
endfunction()
