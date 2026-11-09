# The sources of FLARE's structures. ERF includes this file for each directory in
# ERF_STRUCTURE_MODEL_DIRS (Structures/ERFStructures.cmake adds this one) and calls
# erf_structure_model_add with ERF's OBJECT library; ERF_STRUCTURE_MODEL_DIR holds this
# directory during the call. The conductor model (Conductors/ERF_ConductorModel.cpp)
# registers on ERF's structure coupling as erf.structures.model = conductors.
function(erf_structure_model_add target)
  if(NOT ERF_ENABLE_STRUCTURES)
    message(FATAL_ERROR "FLARE's Structures/ directory is compiled in with ERF_ENABLE_STRUCTURES=ON, "
                        "not through ERF_STRUCTURE_MODEL_DIRS (it needs the MoorDyn set-up of Structures/ERFStructures.cmake)")
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
    ${ERF_STRUCTURE_MODEL_DIR}/Conductors/ERF_GroundContact.cpp
    ${ERF_STRUCTURE_MODEL_DIR}/Conductors/ERF_Conductors.cpp
    ${ERF_STRUCTURE_MODEL_DIR}/Conductors/ERF_ConductorModel.cpp
  )
  target_include_directories(${target} PUBLIC
    $<BUILD_INTERFACE:${ERF_STRUCTURE_MODEL_DIR}/MoorDyn>
    $<BUILD_INTERFACE:${ERF_STRUCTURE_MODEL_DIR}/Towers>
    $<BUILD_INTERFACE:${ERF_STRUCTURE_MODEL_DIR}/Conductors>
    $<INSTALL_INTERFACE:include>)
  # the headers go with ERF's own (CMakeLists.txt installs Source/*.H to include)
  if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    file(GLOB _erf_structures_headers ${ERF_STRUCTURE_MODEL_DIR}/MoorDyn/*.H ${ERF_STRUCTURE_MODEL_DIR}/Towers/*.H
                                       ${ERF_STRUCTURE_MODEL_DIR}/Conductors/*.H)
    install(FILES ${_erf_structures_headers} DESTINATION include)
  endif()
  # read by ERF_MoorDynSystem.cpp only
  target_compile_definitions(${target} PRIVATE ERF_MOORDYN_VERSION="${ERF_MOORDYN_VERSION}")
  if(ERF_MOORDYN_USE_STUB)
    target_compile_definitions(${target} PRIVATE ERF_MOORDYN_USE_STUB)
    target_link_libraries(${target} PUBLIC erf_moordyn_stub)
  else()
    target_link_libraries(${target} PUBLIC MoorDyn::moordyn)
  endif()
endfunction()
