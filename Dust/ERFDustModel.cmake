# The sources of FLARE's dust layer. FLARE's CMakeLists.txt (flare_add_models) includes
# this file with ERF_ENABLE_DUST and calls erf_dust_model_add with the target to compile
# the layer into; ERF_DUST_MODEL_DIR holds this directory during the call. The layer
# registers on ERF's dust coupling as erf.dust.model = flare (ERF_FlareDustModel.cpp),
# the only model in flare and so the one a deck gets when it names none.
function(erf_dust_model_add target)
  if(NOT ERF_ENABLE_FIRE)
    # the dust layer applies the native fire engine's terrain wind corrections
    target_sources(${target} PRIVATE ${PROJECT_SOURCE_DIR}/FireEngine/ERF_FireWindAdjust.cpp)
    target_include_directories(${target} PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/FireEngine>)
  endif()
  target_sources(${target} PRIVATE
    ${ERF_DUST_MODEL_DIR}/ERF_FlareDustModel.cpp
    ${ERF_DUST_MODEL_DIR}/ERF_DustGrid.cpp
    ${ERF_DUST_MODEL_DIR}/ERF_DustPrerequisites.cpp
    ${ERF_DUST_MODEL_DIR}/ERF_DustLayer.cpp
    ${ERF_DUST_MODEL_DIR}/ERF_DustTerrainSlope.cpp
    ${ERF_DUST_MODEL_DIR}/ERF_DustSurfaceReader.cpp
    ${ERF_DUST_MODEL_DIR}/ERF_PhreeqcReader.cpp
    ${ERF_DUST_MODEL_DIR}/ERF_DustEmission.cpp
    ${ERF_DUST_MODEL_DIR}/ERF_DustWindExtract.cpp
    ${ERF_DUST_MODEL_DIR}/ERF_DustAtmCoupling.cpp
    ${ERF_DUST_MODEL_DIR}/ERF_DustPlotfile.cpp
    ${ERF_DUST_MODEL_DIR}/ERF_DustPC.cpp
  )
  if(ERF_ENABLE_FIRE)
    target_sources(${target} PRIVATE
      ${ERF_DUST_MODEL_DIR}/ERF_DustFireLofting.cpp
      ${ERF_DUST_MODEL_DIR}/ERF_FireDustCoupling.cpp
    )
  endif()
  # PUBLIC: the dust unit tests include the layer's headers
  target_include_directories(${target} PUBLIC $<BUILD_INTERFACE:${ERF_DUST_MODEL_DIR}>)
endfunction()
