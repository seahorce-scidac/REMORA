# Unit tests: plain executables that return nonzero on failure, no framework.
# Each gets the "unit" ctest label, so `ctest -L unit` runs them alone and
# `ctest -LE unit` skips them.
set(remora_unit_test_exe_name remora_unit_tests)

add_executable(${remora_unit_test_exe_name}
  ${CMAKE_CURRENT_SOURCE_DIR}/REMORA_DateClock_test.cpp
)

# The lat-lon forcing arithmetic, in every build: the header it tests needs no NetCDF.
set(remora_latlon_test_exe_name remora_latlon_regrid_tests)
add_executable(${remora_latlon_test_exe_name}
  ${CMAKE_CURRENT_SOURCE_DIR}/REMORA_LatLonRegrid_test.cpp
)

foreach (exe ${remora_unit_test_exe_name} ${remora_latlon_test_exe_name})
  target_include_directories(${exe} PRIVATE
    ${PROJECT_SOURCE_DIR}/Source/Utils
    ${PROJECT_SOURCE_DIR}/Source/IO
  )

  # AMReX supplies the inline and GPU qualifiers and the defines that go with
  # them.
  target_link_libraries(${exe} PRIVATE ${REMORA_AMREX_TARGET})

  # AMREX_GPU_HOST_DEVICE expands to __host__ __device__ in a GPU build, so this
  # has to go through the device compiler, as Source/ does.
  if (REMORA_ENABLE_CUDA)
    setup_target_for_cuda_compilation(${exe})
  endif ()
endforeach ()

add_test(NAME DateClock COMMAND ${remora_unit_test_exe_name})
set_tests_properties(DateClock PROPERTIES
  TIMEOUT 60
  LABELS "unit"
)

add_test(NAME LatLonRegrid COMMAND ${remora_latlon_test_exe_name})
set_tests_properties(LatLonRegrid PROPERTIES
  TIMEOUT 60
  LABELS "unit"
)

# The NetCDF readers against every file format the linked PnetCDF can create,
# NetCDF-4 included when PnetCDF has its NetCDF-4 driver. Two ranks, because the
# per-rank reader is collective and is tested with a rank that has nothing to read.
if (REMORA_ENABLE_PNETCDF)
  set(remora_ncformats_test_exe_name remora_ncformats_tests)
  add_executable(${remora_ncformats_test_exe_name}
    ${CMAKE_CURRENT_SOURCE_DIR}/REMORA_NCFormats_test.cpp
    ${PROJECT_SOURCE_DIR}/Source/IO/REMORA_NCInterface.cpp
    ${PROJECT_SOURCE_DIR}/Source/IO/REMORA_NCFile.cpp
  )
  target_include_directories(${remora_ncformats_test_exe_name} PRIVATE
    ${PROJECT_SOURCE_DIR}/Source/Utils
    ${PROJECT_SOURCE_DIR}/Source/IO
    ${PNETCDF_INCLUDE_DIRS}
  )
  target_compile_definitions(${remora_ncformats_test_exe_name} PRIVATE REMORA_USE_NETCDF)
  target_link_libraries(${remora_ncformats_test_exe_name} PRIVATE
    ${REMORA_AMREX_TARGET} ${PNETCDF_LINK_LIBRARIES})
  if (REMORA_ENABLE_CUDA)
    setup_target_for_cuda_compilation(${remora_ncformats_test_exe_name})
  endif ()

  set(remora_ncformats_fixture
      ${PROJECT_SOURCE_DIR}/Tests/test_files/NCFormats/nc_format_test_deflate.nc)
  if (REMORA_ENABLE_MPI)
    add_test(NAME NCFormats
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
              $<TARGET_FILE:${remora_ncformats_test_exe_name}> ${remora_ncformats_fixture})
  else ()
    add_test(NAME NCFormats
      COMMAND ${remora_ncformats_test_exe_name} ${remora_ncformats_fixture})
  endif ()
  set_tests_properties(NCFormats PROPERTIES
    TIMEOUT 60
    LABELS "unit"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
endif ()
