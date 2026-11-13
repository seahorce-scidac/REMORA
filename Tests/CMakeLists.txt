if(WIN32)
    set(FCOMPARE_EXE ${CMAKE_BINARY_DIR}/Submodules/AMReX/Tools/Plotfile/*/amrex_fcompare.exe CACHE INTERNAL "Path to fcompare executable for regression tests")
    set(FEXTREMA_EXE ${CMAKE_BINARY_DIR}/Submodules/AMReX/Tools/Plotfile/*/amrex_fextrema.exe CACHE INTERNAL "Path to fextrema executable for regression tests")
else()
    set(FCOMPARE_EXE ${CMAKE_BINARY_DIR}/Submodules/AMReX/Tools/Plotfile/amrex_fcompare CACHE INTERNAL "Path to fcompare executable for regression tests")
    set(FEXTREMA_EXE ${CMAKE_BINARY_DIR}/Submodules/AMReX/Tools/Plotfile/amrex_fextrema CACHE INTERNAL "Path to fextrema executable for regression tests")
endif()
# Compares a run with dropped boxes (remora.drop_dry_boxes) against a full-layout run, which
# amrex_fcompare cannot do: their BoxArrays cover different cells.
add_executable(remora_subset_compare ${CMAKE_CURRENT_SOURCE_DIR}/tools/remora_subset_compare.cpp)
target_link_libraries(remora_subset_compare PRIVATE ${REMORA_AMREX_TARGET})
if (REMORA_ENABLE_CUDA)
  setup_target_for_cuda_compilation(remora_subset_compare)
endif ()
# Multi-config generators (MSVC) put the executable in a per-configuration subdirectory; the
# tests run it through sh, which expands the glob, as for FCOMPARE_EXE above
if(WIN32)
    set(SUBSET_COMPARE_EXE ${CMAKE_CURRENT_BINARY_DIR}/*/remora_subset_compare.exe)
else()
    set(SUBSET_COMPARE_EXE ${CMAKE_CURRENT_BINARY_DIR}/remora_subset_compare${CMAKE_EXECUTABLE_SUFFIX})
endif()

if(REMORA_BUILD_EXECUTABLES)
    add_executable(remora_bottom_exchange_api
        tools/remora_bottom_exchange_api.cpp
        ${PROJECT_SOURCE_DIR}/Exec/REMORA_Prob.cpp)
    target_include_directories(remora_bottom_exchange_api PRIVATE ${PROJECT_SOURCE_DIR}/Exec)
    target_link_libraries(remora_bottom_exchange_api PRIVATE remora_api)
    add_executable(remora_bottom_exchange_check tools/remora_bottom_exchange_check.cpp)
    target_link_libraries(remora_bottom_exchange_check PRIVATE ${REMORA_AMREX_TARGET})
    if(REMORA_ENABLE_PNETCDF)
        target_compile_definitions(remora_bottom_exchange_check PRIVATE REMORA_USE_NETCDF)
        target_include_directories(remora_bottom_exchange_check PRIVATE ${PNETCDF_INCLUDE_DIRS})
        target_link_libraries(remora_bottom_exchange_check PRIVATE ${PNETCDF_LINK_LIBRARIES})
    endif()
    if(REMORA_ENABLE_CUDA)
        setup_target_for_cuda_compilation(remora_bottom_exchange_api)
        setup_target_for_cuda_compilation(remora_bottom_exchange_check)
    endif()
endif()

add_subdirectory(Unit)

include(${CMAKE_CURRENT_SOURCE_DIR}/CTestList.cmake)
