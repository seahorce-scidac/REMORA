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
