function(target_link_libraries_system target visibility)
  set(libs ${ARGN})
  foreach(lib ${libs})
    get_target_property(lib_include_dirs ${lib} INTERFACE_INCLUDE_DIRECTORIES)
    target_include_directories(${target} SYSTEM ${visibility} ${lib_include_dirs})
    target_link_libraries(${target} ${visibility} ${lib})
  endforeach(lib)
endfunction(target_link_libraries_system)

function(build_remora_lib remora_lib_name)

  set(SRC_DIR ${PROJECT_SOURCE_DIR}/Source)
  set(BIN_DIR ${CMAKE_BINARY_DIR}/Source/${remora_lib_name})

  include(${PROJECT_SOURCE_DIR}/CMake/SetREMORACompileFlags.cmake)
  set_remora_compile_flags(${remora_lib_name})

  if(REMORA_ENABLE_PARTICLES)
    target_sources(${remora_lib_name} PRIVATE
                   ${SRC_DIR}/Particles/REMORA_PC_Evolve.cpp
                   ${SRC_DIR}/Particles/REMORA_PC_Init.cpp
                   ${SRC_DIR}/Particles/REMORA_PC_Utils.cpp
                   ${SRC_DIR}/Particles/REMORA_Tracers.cpp)
    target_include_directories(${remora_lib_name} PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/Source/Particles>)
    target_compile_definitions(${remora_lib_name} PUBLIC REMORA_USE_PARTICLES)
  endif()

  if(REMORA_ENABLE_PNETCDF)
    target_sources(${remora_lib_name} PRIVATE
                   ${SRC_DIR}/IO/REMORA_NCInterface.H
                   ${SRC_DIR}/IO/REMORA_NCPlotFile.H
                   ${SRC_DIR}/IO/REMORA_NCFile.H
                   ${SRC_DIR}/IO/REMORA_NCInterface.cpp
                   ${SRC_DIR}/IO/REMORA_NCPlotFile.cpp
                   ${SRC_DIR}/IO/REMORA_NCFile.cpp
                   ${SRC_DIR}/IO/REMORA_NCTimeSeries.cpp
                   ${SRC_DIR}/IO/REMORA_NCTimeSeriesRiver.cpp
                   ${SRC_DIR}/IO/REMORA_NCTimeSeriesBoundary.cpp
                   ${SRC_DIR}/IO/REMORA_ReadFromInitNetcdf.cpp
                   ${SRC_DIR}/BoundaryConditions/REMORA_BoundaryConditions_netcdf.cpp
                   ${SRC_DIR}/Initialization/REMORA_init_from_netcdf.cpp)
    target_compile_definitions(${remora_lib_name} PUBLIC REMORA_USE_NETCDF)
  endif()

  if(REMORA_ENABLE_MOAB)
    target_sources(${remora_lib_name} PRIVATE
                   ${SRC_DIR}/REMORA_MOAB.cpp
                   ${SRC_DIR}/REMORA_MOAB.H)
    target_compile_definitions(${remora_lib_name} PUBLIC REMORA_USE_MOAB)
  endif()

  if(REMORA_ENABLE_FUNWAVE_FORT)
  target_sources(${remora_lib_name}
     PRIVATE
#      ${SRC_DIR}/FunWave/REMORA_funwave.F90
       )
  target_compile_definitions(${remora_lib_name} PUBLIC REMORA_USE_FUNWAVE_FORT)
  endif()

  if((NOT REMORA_BUILD_EXECUTABLES) OR REMORA_BUILD_LIBRARY_ONLY)
    # In library-only superbuild mode, archive extraction + weak amrex_probinit
    # requires a forced reference path (see REMORA.cpp/REMORA_Prob.cpp link anchor).
    # Avoid cross-application symbol collision when REMORA and ERF are linked
    # into one parent executable by forcing explicit REMORA-prefixed names.
    target_compile_definitions(${remora_lib_name} PRIVATE
                   ERF_REMORA_FORCE_PROBINIT_LINK=1
                   amrex_probinit=remora_probinit)
    target_compile_definitions(${remora_lib_name} PRIVATE
                   Problem=REMORAProblem
                   ProblemBase=REMORAProblemBase
                   SolverChoice=REMORASolverChoice)
    target_sources(${remora_lib_name} PRIVATE
                   ${PROJECT_SOURCE_DIR}/Exec/REMORA_Prob.cpp)
  endif()

  # Coupling source is present only on coupling branches.
  # Build/link branches should compile without requiring this file.
  if(EXISTS "${SRC_DIR}/REMORA_Coupling.cpp")
    target_sources(${remora_lib_name} PRIVATE
                   ${SRC_DIR}/REMORA_Coupling.cpp)
  endif()

  target_sources(${remora_lib_name}
     PRIVATE
       ${SRC_DIR}/REMORA_Derive.cpp
       ${SRC_DIR}/REMORA.cpp
       ${SRC_DIR}/REMORA_SumIQ.cpp
       ${SRC_DIR}/REMORA_Tagging.cpp
       ${SRC_DIR}/BoundaryConditions/REMORA_BoundaryConditions_cons.cpp
       ${SRC_DIR}/BoundaryConditions/REMORA_BoundaryConditions_xvel.cpp
       ${SRC_DIR}/BoundaryConditions/REMORA_BoundaryConditions_yvel.cpp
       ${SRC_DIR}/BoundaryConditions/REMORA_BoundaryConditions_zvel.cpp
       ${SRC_DIR}/BoundaryConditions/REMORA_FillPatch.cpp
       ${SRC_DIR}/BoundaryConditions/REMORA_FillPatcher.cpp
       ${SRC_DIR}/BoundaryConditions/REMORA_PhysBCFunct.cpp
       ${SRC_DIR}/Initialization/REMORA_init.cpp
       ${SRC_DIR}/Initialization/REMORA_init1d.cpp
       ${SRC_DIR}/Initialization/REMORA_init_bcs.cpp
       ${SRC_DIR}/Initialization/REMORA_make_new_level.cpp
       ${SRC_DIR}/IO/REMORA_Checkpoint.cpp
       ${SRC_DIR}/IO/REMORA_Plotfile.cpp
       ${SRC_DIR}/IO/REMORA_SetPlotVars.cpp
       ${SRC_DIR}/IO/REMORA_writeJobInfo.cpp
       ${SRC_DIR}/IO/REMORA_console_io.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_Advance.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_advance_2d.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_advance_2d_onestep.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_advance_3d.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_advance_3d_ml.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_apply_clim_nudg.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_bulk_flux.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_biology.cpp
      ${SRC_DIR}/TimeIntegration/Biology/REMORA_bio_fennel.cpp
      ${SRC_DIR}/TimeIntegration/Biology/REMORA_bio_placeholder.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_setup_step.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_rho_eos.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_gls.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_prsgrd.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_uv3dmix.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_t3dmix.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_coriolis.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_prestep.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_prestep_t_advection.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_prestep_diffusion.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_rhs_t_3d.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_rhs_uv_3d.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_rhs_uv_2d.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_scale_rhs_vars.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_vert_visc_3d.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_update_massflux_3d.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_vert_mean_3d.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_ComputeTimestep.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_TimeStep.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_TimeStepML.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_set_weights.cpp
  )

  # Example executables rely on Source/main.cpp for the process entry point
  # and for globally shared CLI metadata (e.g., inputs_name).
  if(REMORA_BUILD_EXECUTABLES)
    target_sources(${remora_lib_name}
       PRIVATE
         ${SRC_DIR}/main.cpp
    )
  endif()

  include(AMReXBuildInfo)
  generate_buildinfo(${remora_lib_name} ${PROJECT_SOURCE_DIR})
  if(AMREX_C_SCRIPTS_DIR)
    target_include_directories(${remora_lib_name} PUBLIC $<BUILD_INTERFACE:${AMREX_C_SCRIPTS_DIR}>)
  elseif(REMORA_USE_INTERNAL_AMREX)
    target_include_directories(${remora_lib_name} PUBLIC $<BUILD_INTERFACE:${AMREX_SUBMOD_LOCATION}/Tools/C_scripts>)
  endif()

  if(REMORA_ENABLE_PNETCDF)
    if(PNETCDF_FOUND)
      #Link our executable to the PNETCDF libraries, etc
      target_link_libraries(${remora_lib_name} PUBLIC ${PNETCDF_LINK_LIBRARIES})
      target_include_directories(${remora_lib_name} PUBLIC ${PNETCDF_INCLUDE_DIRS})
    endif()
  endif()

  if(REMORA_ENABLE_MPI)
    target_link_libraries(${remora_lib_name} PUBLIC $<$<BOOL:${MPI_CXX_FOUND}>:MPI::MPI_CXX>)
  endif()

  if(REMORA_ENABLE_MOAB)
    if(MOAB_FOUND)
      #Link our executable to the MOAB libraries, etc
      target_link_libraries(${remora_lib_name} PUBLIC ${MOAB_LIBRARIES})
      target_include_directories(${remora_lib_name} PUBLIC ${MOAB_INCLUDE_DIRS})
    endif()
  endif()

  #REMORA include directories
  target_include_directories(${remora_lib_name} PUBLIC  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/Source>)
  target_include_directories(${remora_lib_name} PUBLIC  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/Source/BoundaryConditions>)
  target_include_directories(${remora_lib_name} PUBLIC  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/Source/Initialization>)
  target_include_directories(${remora_lib_name} PUBLIC  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/Source/Utils>)
  target_include_directories(${remora_lib_name} PUBLIC  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/Source/TimeIntegration>)
  target_include_directories(${remora_lib_name} PUBLIC  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/Source/IO>)
  target_include_directories(${remora_lib_name} PUBLIC  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>)

  #Link to amrex library
  target_link_libraries_system(${remora_lib_name} PUBLIC ${REMORA_AMREX_TARGET})
  if(REMORA_ENABLE_CUDA)
    set(pctargets "${remora_lib_name}")
    foreach(tgt IN LISTS pctargets)
      get_target_property(REMORA_SOURCES ${tgt} SOURCES)
      list(FILTER REMORA_SOURCES INCLUDE REGEX "\\.cpp")
      set_source_files_properties(${REMORA_SOURCES} PROPERTIES LANGUAGE CUDA)
      message(STATUS "setting cuda for ${REMORA_SOURCES}")
    endforeach()
    set_target_properties(
    ${remora_lib_name} PROPERTIES
    LANGUAGE CUDA
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_RESOLVE_DEVICE_SYMBOLS ON)
  endif()

  #Define what we want to be installed during a make install
  install(TARGETS ${remora_lib_name}
          RUNTIME DESTINATION bin
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib)

endfunction(build_remora_lib)

function(build_remora_exe remora_exe_name)

  set(SRC_DIR ${PROJECT_SOURCE_DIR}/Source)

  target_link_libraries(${remora_exe_name}  PUBLIC ${remora_lib_name})
  include(${PROJECT_SOURCE_DIR}/CMake/SetREMORACompileFlags.cmake)
  set_remora_compile_flags(${remora_exe_name})

  if(REMORA_ENABLE_CUDA)
    set(pctargets "${remora_exe_name}")
    foreach(tgt IN LISTS pctargets)
      get_target_property(REMORA_SOURCES ${tgt} SOURCES)
      list(FILTER REMORA_SOURCES INCLUDE REGEX "\\.cpp")
      set_source_files_properties(${REMORA_SOURCES} PROPERTIES LANGUAGE CUDA)
      message(STATUS "setting cuda for ${REMORA_SOURCES}")
    endforeach()
    set_target_properties(
    ${remora_exe_name} PROPERTIES
    LANGUAGE CUDA
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_RESOLVE_DEVICE_SYMBOLS ON)
  endif()

  install(TARGETS ${remora_exe_name}
          RUNTIME DESTINATION bin
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib)

endfunction()
