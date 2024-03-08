function(target_link_libraries_system target visibility)
  set(libs ${ARGN})
  foreach(lib ${libs})
    get_target_property(lib_include_dirs ${lib} INTERFACE_INCLUDE_DIRECTORIES)
    target_include_directories(${target} SYSTEM ${visibility} ${lib_include_dirs})
    target_link_libraries(${target} ${visibility} ${lib})
  endforeach(lib)
endfunction(target_link_libraries_system)

function(build_remora_lib remora_lib_name)

  set(SRC_DIR ${CMAKE_SOURCE_DIR}/Source)
  set(BIN_DIR ${CMAKE_BINARY_DIR}/Source/${remora_lib_name})

  include(${CMAKE_SOURCE_DIR}/CMake/SetREMORACompileFlags.cmake)
  set_remora_compile_flags(${remora_lib_name})

  set(REMORA_EOS_DIR "${CMAKE_SOURCE_DIR}/Source")
  target_sources(${remora_lib_name} PRIVATE
                 ${REMORA_EOS_DIR}/EOS.H)
  target_include_directories(${remora_lib_name} SYSTEM PUBLIC ${REMORA_EOS_DIR})

  if(REMORA_ENABLE_PARTICLES)
    target_sources(${remora_lib_name} PRIVATE
                   ${SRC_DIR}/Particles/TracerPC.cpp)
    target_include_directories(${remora_lib_name} PUBLIC ${SRC_DIR}/Particles)
    target_compile_definitions(${remora_lib_name} PUBLIC REMORA_USE_PARTICLES)
  endif()

  if(REMORA_ENABLE_NETCDF)
    target_sources(${remora_lib_name} PRIVATE
                   ${SRC_DIR}/IO/NCInterface.H
                   ${SRC_DIR}/IO/NCPlotFile.H
                   ${SRC_DIR}/IO/NCFile.H
                   ${SRC_DIR}/IO/NCInterface.cpp
                   ${SRC_DIR}/IO/NCPlotFile.cpp
                   ${SRC_DIR}/IO/NCFile.cpp
                   ${SRC_DIR}/IO/ReadFromInitNetcdf.cpp
                   ${SRC_DIR}/IO/ReadFromBdryNetcdf.cpp
                   ${SRC_DIR}/BoundaryConditions/BoundaryConditions_netcdf.cpp
                   ${SRC_DIR}/Initialization/REMORA_init_from_netcdf.cpp)
    target_compile_definitions(${remora_lib_name} PUBLIC REMORA_USE_NETCDF)
  endif()

  if(REMORA_ENABLE_HDF5)
    target_compile_definitions(${remora_lib_name} PUBLIC REMORA_USE_HDF5)
  endif()

  target_sources(${remora_lib_name}
     PRIVATE
       ${SRC_DIR}/Derive.cpp
       ${SRC_DIR}/REMORA.cpp
       ${SRC_DIR}/REMORA_SumIQ.cpp
       ${SRC_DIR}/REMORA_Tagging.cpp
       ${SRC_DIR}/BoundaryConditions/BoundaryConditions_cons.cpp
       ${SRC_DIR}/BoundaryConditions/BoundaryConditions_xvel.cpp
       ${SRC_DIR}/BoundaryConditions/BoundaryConditions_yvel.cpp
       ${SRC_DIR}/BoundaryConditions/BoundaryConditions_zvel.cpp
       ${SRC_DIR}/BoundaryConditions/REMORA_FillPatch.cpp
       ${SRC_DIR}/BoundaryConditions/REMORA_PhysBCFunct.cpp
       ${SRC_DIR}/Initialization/REMORA_init.cpp
       ${SRC_DIR}/Initialization/REMORA_init1d.cpp
       ${SRC_DIR}/Initialization/REMORA_init_bcs.cpp
       ${SRC_DIR}/Initialization/REMORA_make_new_level.cpp
       ${SRC_DIR}/IO/Checkpoint.cpp
       ${SRC_DIR}/IO/Plotfile.cpp
       ${SRC_DIR}/IO/writeJobInfo.cpp
       ${SRC_DIR}/IO/console_io.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_Advance.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_advance_2d.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_advance_2d_onestep.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_advance_3d.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_advance_3d_ml.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_setup_step.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_rho_eos.cpp
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
       ${SRC_DIR}/TimeIntegration/REMORA_vert_visc_3d.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_update_massflux_3d.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_vert_mean_3d.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_ComputeTimestep.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_TimeStep.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_TimeStepML.cpp
       ${SRC_DIR}/TimeIntegration/REMORA_set_weights.cpp
  )

  if(NOT "${remora_exe_name}" STREQUAL "remora_unit_tests")
    target_sources(${remora_lib_name}
       PRIVATE
         ${SRC_DIR}/main.cpp
    )
  endif()

  include(AMReXBuildInfo)
  generate_buildinfo(${remora_lib_name} ${CMAKE_SOURCE_DIR})
  target_include_directories(${remora_lib_name} PUBLIC ${AMREX_SUBMOD_LOCATION}/Tools/C_scripts)

  if(REMORA_ENABLE_NETCDF)
    if(NETCDF_FOUND)
      #Link our executable to the NETCDF libraries, etc
      target_link_libraries(${remora_lib_name} PUBLIC ${NETCDF_LINK_LIBRARIES})
      target_include_directories(${remora_lib_name} PUBLIC ${NETCDF_INCLUDE_DIRS})
    endif()
  endif()

  if(REMORA_ENABLE_MPI)
    target_link_libraries(${remora_lib_name} PUBLIC $<$<BOOL:${MPI_CXX_FOUND}>:MPI::MPI_CXX>)
  endif()

  #REMORA include directories
  target_include_directories(${remora_lib_name} PUBLIC ${SRC_DIR})
  target_include_directories(${remora_lib_name} PUBLIC ${SRC_DIR}/BoundaryConditions)
  target_include_directories(${remora_lib_name} PUBLIC ${SRC_DIR}/Initialization)
  target_include_directories(${remora_lib_name} PUBLIC ${SRC_DIR}/Utils)
  target_include_directories(${remora_lib_name} PUBLIC ${SRC_DIR}/TimeIntegration)
  target_include_directories(${remora_lib_name} PUBLIC ${SRC_DIR}/IO)
  target_include_directories(${remora_lib_name} PUBLIC ${CMAKE_BINARY_DIR})

  #Link to amrex library
  target_link_libraries_system(${remora_lib_name} PUBLIC amrex)
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

  set(SRC_DIR ${CMAKE_SOURCE_DIR}/Source)

  target_link_libraries(${remora_exe_name}  PUBLIC ${remora_lib_name})
  include(${CMAKE_SOURCE_DIR}/CMake/SetREMORACompileFlags.cmake)
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
