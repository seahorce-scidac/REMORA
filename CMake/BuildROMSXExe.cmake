function(target_link_libraries_system target visibility)
  set(libs ${ARGN})
  foreach(lib ${libs})
    get_target_property(lib_include_dirs ${lib} INTERFACE_INCLUDE_DIRECTORIES)
    target_include_directories(${target} SYSTEM ${visibility} ${lib_include_dirs})
    target_link_libraries(${target} ${visibility} ${lib})
  endforeach(lib)
endfunction(target_link_libraries_system)

function(build_romsx_lib romsx_lib_name)

  set(SRC_DIR ${CMAKE_SOURCE_DIR}/Source)
  set(BIN_DIR ${CMAKE_BINARY_DIR}/Source/${romsx_lib_name})

  include(${CMAKE_SOURCE_DIR}/CMake/SetROMSXCompileFlags.cmake)
  set_romsx_compile_flags(${romsx_lib_name})

  set(ROMSX_EOS_DIR "${CMAKE_SOURCE_DIR}/Source")
  target_sources(${romsx_lib_name} PRIVATE
                 ${ROMSX_EOS_DIR}/EOS.H)
  target_include_directories(${romsx_lib_name} SYSTEM PUBLIC ${ROMSX_EOS_DIR})

  if(ROMSX_ENABLE_NETCDF)
    target_sources(${romsx_lib_name} PRIVATE
                   ${SRC_DIR}/IO/NCIntromsxace.H
                   ${SRC_DIR}/IO/NCWpsFile.H
                   ${SRC_DIR}/IO/NCPlotFile.H
                   ${SRC_DIR}/IO/NCBuildFABs.cpp
                   ${SRC_DIR}/IO/NCIntromsxace.cpp
                   ${SRC_DIR}/IO/NCPlotFile.cpp
                   ${SRC_DIR}/IO/NCCheckpoint.cpp
                   ${SRC_DIR}/IO/NCMultiFabFile.cpp)
    target_compile_definitions(${romsx_lib_name} PUBLIC ROMSX_USE_NETCDF)
  endif()

  if(ROMSX_ENABLE_HDF5)
    target_compile_definitions(${romsx_lib_name} PUBLIC ROMSX_USE_HDF5)
  endif()

  target_sources(${romsx_lib_name}
     PRIVATE
       ${SRC_DIR}/DataStruct.H
       ${SRC_DIR}/ROMSX_Constants.H
       ${SRC_DIR}/Derive.H
       ${SRC_DIR}/Derive.cpp
       ${SRC_DIR}/IndexDefines.H
       ${SRC_DIR}/prob_common.H
       ${SRC_DIR}/ROMSX.H
       ${SRC_DIR}/ROMSX.cpp
       ${SRC_DIR}/ROMSX_init.cpp
       ${SRC_DIR}/ROMSX_init1d.cpp
       ${SRC_DIR}/ROMSX_SumIQ.cpp
       ${SRC_DIR}/ROMSX_Tagging.cpp
       ${SRC_DIR}/BoundaryConditions/BoundaryConditions_cons.cpp
       ${SRC_DIR}/BoundaryConditions/BoundaryConditions_xvel.cpp
       ${SRC_DIR}/BoundaryConditions/BoundaryConditions_yvel.cpp
       ${SRC_DIR}/BoundaryConditions/BoundaryConditions_zvel.cpp
       ${SRC_DIR}/BoundaryConditions/ROMSX_FillPatch.cpp
       ${SRC_DIR}/BoundaryConditions/ROMSX_PhysBCFunct.cpp
       ${SRC_DIR}/BoundaryConditions/PlaneAverage.H
       ${SRC_DIR}/BoundaryConditions/VelPlaneAverage.H
       ${SRC_DIR}/BoundaryConditions/DirectionSelector.H
       ${SRC_DIR}/IO/Checkpoint.cpp
       ${SRC_DIR}/IO/Plotfile.cpp
       ${SRC_DIR}/IO/writeJobInfo.cpp
       ${SRC_DIR}/Utils/ROMSX_Math.H
       ${SRC_DIR}/Utils/Interpolation.H
       ${SRC_DIR}/Utils/MomentumToVelocity.cpp
       ${SRC_DIR}/Utils/Utils.H
       ${SRC_DIR}/Utils/TerrainMetrics.H
       ${SRC_DIR}/Utils/TerrainMetrics.cpp
       ${SRC_DIR}/Utils/VelocityToMomentum.cpp
       ${SRC_DIR}/TimeIntegration/ROMSX_Advance.cpp
       ${SRC_DIR}/TimeIntegration/ROMSX_advance_2d.cpp
       ${SRC_DIR}/TimeIntegration/ROMSX_advance_3d.cpp
       ${SRC_DIR}/TimeIntegration/ROMSX_coriolis.cpp
       ${SRC_DIR}/TimeIntegration/ROMSX_prestep_t_3d.cpp
       ${SRC_DIR}/TimeIntegration/ROMSX_prestep_uv_3d.cpp
       ${SRC_DIR}/TimeIntegration/ROMSX_rhs_t_3d.cpp
       ${SRC_DIR}/TimeIntegration/ROMSX_rhs_uv_3d.cpp
       ${SRC_DIR}/TimeIntegration/ROMSX_update_vel_3d.cpp
       ${SRC_DIR}/TimeIntegration/ROMSX_vert_visc_3d.cpp
       ${SRC_DIR}/TimeIntegration/ROMSX_set_massflux_3d.cpp
       ${SRC_DIR}/TimeIntegration/ROMSX_update_massflux_3d.cpp
       ${SRC_DIR}/TimeIntegration/ROMSX_vert_mean_3d.cpp
       ${SRC_DIR}/TimeIntegration/ROMSX_ComputeTimestep.cpp
       ${SRC_DIR}/TimeIntegration/ROMSX_TimeStep.cpp
  )

  if(NOT "${romsx_exe_name}" STREQUAL "romsx_unit_tests")
    target_sources(${romsx_lib_name}
       PRIVATE
         ${SRC_DIR}/main.cpp
    )
  endif()

  include(AMReXBuildInfo)
  generate_buildinfo(${romsx_lib_name} ${CMAKE_SOURCE_DIR})
  target_include_directories(${romsx_lib_name} PUBLIC ${AMREX_SUBMOD_LOCATION}/Tools/C_scripts)

  if(ROMSX_ENABLE_NETCDF)
    if(NETCDF_FOUND)
      #Link our executable to the NETCDF libraries, etc
      target_link_libraries(${romsx_lib_name} PUBLIC ${NETCDF_LINK_LIBRARIES})
      target_include_directories(${romsx_lib_name} PUBLIC ${NETCDF_INCLUDE_DIRS})
    endif()
  endif()

  if(ROMSX_ENABLE_MPI)
    target_link_libraries(${romsx_lib_name} PUBLIC $<$<BOOL:${MPI_CXX_FOUND}>:MPI::MPI_CXX>)
  endif()

  #ROMSX include directories
  target_include_directories(${romsx_lib_name} PUBLIC ${SRC_DIR})
  target_include_directories(${romsx_lib_name} PUBLIC ${SRC_DIR}/BoundaryConditions)
  target_include_directories(${romsx_lib_name} PUBLIC ${SRC_DIR}/Utils)
  target_include_directories(${romsx_lib_name} PUBLIC ${SRC_DIR}/TimeIntegration)
  target_include_directories(${romsx_lib_name} PUBLIC ${SRC_DIR}/IO)
  target_include_directories(${romsx_lib_name} PUBLIC ${CMAKE_BINARY_DIR})

  #Link to amrex library
  target_link_libraries_system(${romsx_lib_name} PUBLIC amrex)
  if(ROMSX_ENABLE_CUDA)
    set(pctargets "${romsx_lib_name}")
    foreach(tgt IN LISTS pctargets)
      get_target_property(ROMSX_SOURCES ${tgt} SOURCES)
      list(FILTER ROMSX_SOURCES INCLUDE REGEX "\\.cpp")
      set_source_files_properties(${ROMSX_SOURCES} PROPERTIES LANGUAGE CUDA)
      message(STATUS "setting cuda for ${ROMSX_SOURCES}")
    endforeach()
    set_target_properties(
    ${romsx_lib_name} PROPERTIES
    LANGUAGE CUDA
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_RESOLVE_DEVICE_SYMBOLS ON)
  endif()

  #Define what we want to be installed during a make install
  install(TARGETS ${romsx_lib_name}
          RUNTIME DESTINATION bin
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib)

endfunction(build_romsx_lib)

function(build_romsx_exe romsx_exe_name)

  set(SRC_DIR ${CMAKE_SOURCE_DIR}/Source)

  target_link_libraries(${romsx_exe_name}  PUBLIC ${romsx_lib_name})
  include(${CMAKE_SOURCE_DIR}/CMake/SetROMSXCompileFlags.cmake)
  set_romsx_compile_flags(${romsx_exe_name})

  target_sources(${romsx_exe_name}
     PRIVATE
       ${SRC_DIR}/ROMSX_init_bcs.cpp

  )
  if(ROMSX_ENABLE_CUDA)
    set(pctargets "${romsx_exe_name}")
    foreach(tgt IN LISTS pctargets)
      get_target_property(ROMSX_SOURCES ${tgt} SOURCES)
      list(FILTER ROMSX_SOURCES INCLUDE REGEX "\\.cpp")
      set_source_files_properties(${ROMSX_SOURCES} PROPERTIES LANGUAGE CUDA)
      message(STATUS "setting cuda for ${ROMSX_SOURCES}")
    endforeach()
    set_target_properties(
    ${romsx_exe_name} PROPERTIES
    LANGUAGE CUDA
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_RESOLVE_DEVICE_SYMBOLS ON)
  endif()

  install(TARGETS ${romsx_exe_name}
          RUNTIME DESTINATION bin
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib)

endfunction()
