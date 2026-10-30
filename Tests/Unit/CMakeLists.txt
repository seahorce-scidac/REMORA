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

# Build the native biology interface independently of the optional Fortran
# oracle. Linking all remora_srclib objects would also bring in main.cpp.
add_executable(remora_biology_tests
  REMORA_Biology_test.cpp
  ${PROJECT_SOURCE_DIR}/Source/Biology/REMORA_Biology.cpp
  ${PROJECT_SOURCE_DIR}/Source/Biology/REMORA_FennelBiology.cpp
)
get_target_property(remora_biology_include_dirs remora_srclib INCLUDE_DIRECTORIES)
target_include_directories(remora_biology_tests PRIVATE ${remora_biology_include_dirs})
target_link_libraries(remora_biology_tests PRIVATE ${REMORA_AMREX_TARGET})
# The parity diagnostics are header-only printf inside the Fennel kernel, so
# compile the test TU the way remora_srclib is compiled. The Fortran bridge
# stays out: it would pull the whole ROMS kernel into the link.
if (REMORA_ENABLE_BIOLOGY_DIAG)
  target_compile_definitions(remora_biology_tests PRIVATE REMORA_USE_BIOLOGY_DIAG)
endif ()
if (REMORA_ENABLE_CUDA)
  setup_target_for_cuda_compilation(remora_biology_tests)
endif ()
add_test(NAME BiologyInterface COMMAND remora_biology_tests)
set_tests_properties(BiologyInterface PROPERTIES TIMEOUT 60 LABELS "unit")
