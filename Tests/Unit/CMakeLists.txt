# Unit tests: plain executables that return nonzero on failure, no framework.
# Each gets the "unit" ctest label, so `ctest -L unit` runs them alone and
# `ctest -LE unit` skips them.
set(remora_unit_test_exe_name remora_unit_tests)

add_executable(${remora_unit_test_exe_name}
  ${CMAKE_CURRENT_SOURCE_DIR}/REMORA_DateClock_test.cpp
)

target_include_directories(${remora_unit_test_exe_name} PRIVATE
  ${PROJECT_SOURCE_DIR}/Source/Utils
)

# AMReX supplies the inline and GPU qualifiers and the defines that go with
# them.
target_link_libraries(${remora_unit_test_exe_name} PRIVATE ${REMORA_AMREX_TARGET})

# AMREX_GPU_HOST_DEVICE expands to __host__ __device__ in a GPU build, so this
# has to go through the device compiler, as Source/ does.
if (REMORA_ENABLE_CUDA)
  setup_target_for_cuda_compilation(${remora_unit_test_exe_name})
endif ()

add_test(NAME DateClock COMMAND ${remora_unit_test_exe_name})
set_tests_properties(DateClock PROPERTIES
  TIMEOUT 60
  LABELS "unit"
)
