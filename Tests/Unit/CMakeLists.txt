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

# Build the native biology interface independently of the optional Fortran
# oracle. Linking all remora_srclib objects would also bring in main.cpp.
add_executable(remora_biology_tests
  REMORA_Biology_test.cpp
  ${PROJECT_SOURCE_DIR}/Source/Biology/REMORA_Biology.cpp
  ${PROJECT_SOURCE_DIR}/Source/Biology/REMORA_FennelBiology.cpp
  ${PROJECT_SOURCE_DIR}/Source/Biology/REMORA_ECBBiology.cpp
)
get_target_property(remora_biology_include_dirs remora_srclib INCLUDE_DIRECTORIES)
target_include_directories(remora_biology_tests PRIVATE ${remora_biology_include_dirs})
target_link_libraries(remora_biology_tests PRIVATE ${REMORA_AMREX_TARGET})
# The parity diagnostics are header-only printf inside the Fennel kernel, so
# compile the test TU the way remora_srclib is compiled. The Fortran bridge
# stays out: it would pull the whole ROMS kernel into the link.
if (REMORA_ENABLE_BIOLOGY_DIAG)
  target_compile_definitions(remora_biology_tests
    PRIVATE REMORA_USE_BIOLOGY_DIAG)
endif ()
if (REMORA_ENABLE_CUDA)
  setup_target_for_cuda_compilation(remora_biology_tests)
endif ()
add_test(NAME BiologyInterface COMMAND remora_biology_tests)
set_tests_properties(BiologyInterface PROPERTIES TIMEOUT 60 LABELS "unit")

if(REMORA_BUILD_EXECUTABLES)
  # Link the archive so Source/main.cpp is not extracted: this driver supplies
  # its own main, inputs_name and analytic problem, and exercises the public API.
  add_executable(remora_bottom_exchange_tests REMORA_BottomExchange_test.cpp)
  target_link_libraries(remora_bottom_exchange_tests PRIVATE remora_api)
  if(REMORA_ENABLE_CUDA)
    setup_target_for_cuda_compilation(remora_bottom_exchange_tests)
  endif()
  function(bottom_exchange_test name)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${name})
    add_test(NAME ${name} COMMAND remora_bottom_exchange_tests
      ${CMAKE_CURRENT_SOURCE_DIR}/bottom_exchange.i ${ARGN})
    set_tests_properties(${name} PROPERTIES TIMEOUT 120 LABELS "bottom_exchange"
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${name})
  endfunction()
  bottom_exchange_test(BottomFreshwater)
  bottom_exchange_test(BottomMatching remora.bottom_salt=35)
  bottom_exchange_test(BottomOutflow remora.bottom_q=-0.0001)
  bottom_exchange_test(BottomImpermeable remora.bottom_exchange_type=none)
  bottom_exchange_test(BottomAnalytic remora.bottom_exchange_type=analytic)
  bottom_exchange_test(BottomCoupler remora.bottom_exchange_type=coupler)
  bottom_exchange_test(BottomGLS remora.vertical_mixing_type=GLS remora.prob.u_0=0
    remora.bottom_exchange_diagnostics=true remora.Scoef=0.00076 test.unstable=true)
  bottom_exchange_test(BottomAMR amr.max_level=1 "amr.ref_ratio_vect=2 2 1"
    remora.refinement_indicators=patch "remora.patch.in_box_lo=200 200"
    "remora.patch.in_box_hi=600 600" remora.patch.max_level=1 remora.do_substep=true
    remora.cf_fill_all_kcomp=0 remora.cf_time_interp_zeta=0 remora.prob.u_0=0)
  bottom_exchange_test(BottomAMRLockstep amr.max_level=1 "amr.ref_ratio_vect=2 2 1"
    remora.refinement_indicators=patch "remora.patch.in_box_lo=200 200"
    "remora.patch.in_box_hi=600 600" remora.patch.max_level=1 remora.do_substep=false
    remora.prob.u_0=0)
  bottom_exchange_test(BottomRegrid amr.max_level=1 "amr.ref_ratio_vect=2 2 1" amr.regrid_int=1
    remora.refinement_indicators=patch "remora.patch.in_box_lo=200 200"
    "remora.patch.in_box_hi=600 600" remora.patch.max_level=1 remora.patch.start_time=2
    remora.cf_fill_all_kcomp=0 remora.cf_time_interp_zeta=0 remora.prob.u_0=0)
  bottom_exchange_test(BottomSpatialAMR remora.bottom_exchange_type=analytic test.spatial=true
    remora.bottom_salt=35 amr.max_level=1 "amr.ref_ratio_vect=2 2 1"
    remora.refinement_indicators=patch "remora.patch.in_box_lo=0 0"
    "remora.patch.in_box_hi=800 800" remora.patch.max_level=1 remora.prob.u_0=0)
  bottom_exchange_test(BottomRestartWrite remora.check_int=2)
  bottom_exchange_test(BottomRestart test.resume=true
    remora.restart=${CMAKE_CURRENT_BINARY_DIR}/BottomRestartWrite/chk00002)
  set_tests_properties(BottomRestart PROPERTIES DEPENDS BottomRestartWrite)
  bottom_exchange_test(BottomReverse remora.bottom_exchange_type=coupler test.reverse=true
    remora.bottom_salt=35 remora.bottom_tracer=2 remora.prob.u_0=0)
  function(bottom_exchange_error name message)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${name})
    add_test(NAME ${name} COMMAND sh ${PROJECT_SOURCE_DIR}/Tests/check_bottom_exchange_error.sh
      "${message}" $<TARGET_FILE:remora_bottom_exchange_tests>
      ${CMAKE_CURRENT_SOURCE_DIR}/bottom_exchange.i ${ARGN})
    set_tests_properties(${name} PROPERTIES TIMEOUT 120 LABELS "bottom_exchange"
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${name})
  endfunction()
  bottom_exchange_error(BottomMissingTracer "explicit finite concentrations" test.omit_salt=true)
  bottom_exchange_error(BottomExpired "forcing is missing or expired"
    remora.bottom_exchange_type=coupler test.valid_end=1)
  bottom_exchange_error(BottomCourant "exceeds bottom_exchange_max_courant" remora.bottom_q=1)
  if(REMORA_ENABLE_PNETCDF)
    foreach(kind 1 2 3)
      bottom_exchange_test(BottomNetCDF${kind} remora.bottom_exchange_type=netcdf
        remora.bottom_exchange_file=bottom.nc test.file_kind=${kind})
    endforeach()
    foreach(mode coupler netcdf)
      set(name BottomMPI_${mode})
      file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${name})
      set(options remora.bottom_exchange_type=${mode} amr.max_grid_size=4)
      if(mode STREQUAL "netcdf")
        list(APPEND options remora.bottom_exchange_file=bottom.nc test.file_kind=3)
      endif()
      add_test(NAME ${name} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
        ${MPIEXEC_PREFLAGS} $<TARGET_FILE:remora_bottom_exchange_tests>
        ${CMAKE_CURRENT_SOURCE_DIR}/bottom_exchange.i ${options} ${MPIEXEC_POSTFLAGS})
      set_tests_properties(${name} PROPERTIES TIMEOUT 120 LABELS "bottom_exchange;mpi"
        PROCESSORS 2 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${name})
    endforeach()
  endif()
endif()
