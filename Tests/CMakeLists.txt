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
set(SUBSET_COMPARE_EXE ${CMAKE_CURRENT_BINARY_DIR}/remora_subset_compare${CMAKE_EXECUTABLE_SUFFIX})

include(${CMAKE_CURRENT_SOURCE_DIR}/CTestList.cmake)

add_subdirectory(Unit)
