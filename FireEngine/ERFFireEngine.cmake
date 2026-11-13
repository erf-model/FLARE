# FLARE's fire engine, for CMake: FLARE's CMakeLists.txt (flare_add_models) includes
# this file and calls erf_fire_engine_add with the target to compile the engine into
# (ERF_FIRE_ENGINE_DIR is this directory during the call). The engine registers on ERF's
# fire coupling as erf.fire.engine = flare.
function(erf_fire_engine_add target)
  target_sources(${target} PRIVATE
    ${ERF_FIRE_ENGINE_DIR}/ERF_NativeFireEngine.cpp
    ${ERF_FIRE_ENGINE_DIR}/ERF_FireEnginePrerequisites.cpp
    ${ERF_FIRE_ENGINE_DIR}/ERF_FireWindAdjust.cpp
    ${ERF_FIRE_ENGINE_DIR}/ERF_FireTerrainReader.cpp
    ${ERF_FIRE_ENGINE_DIR}/ERF_TerrainSlope.cpp
    ${ERF_FIRE_ENGINE_DIR}/ERF_Rothermel.cpp
    ${ERF_FIRE_ENGINE_DIR}/ERF_FireLayer.cpp
    ${ERF_FIRE_ENGINE_DIR}/ERF_FireCheckpoint.cpp
    ${ERF_FIRE_ENGINE_DIR}/ERF_FireSuppression.cpp
    ${ERF_FIRE_ENGINE_DIR}/ERF_Fire.cpp
    ${ERF_FIRE_ENGINE_DIR}/ERF_FirePlotfile.cpp
    ${ERF_FIRE_ENGINE_DIR}/ERF_StructureIgnition.cpp
    ${ERF_FIRE_ENGINE_DIR}/ERF_CustomFuel.cpp
    ${ERF_FIRE_ENGINE_DIR}/ERF_FarsiteWindOutput.cpp
  )
  target_include_directories(${target} PUBLIC $<BUILD_INTERFACE:${ERF_FIRE_ENGINE_DIR}>)
endfunction()
