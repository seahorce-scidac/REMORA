function(set_remora_compile_flags target)
# Logic for handling warnings
if(REMORA_ENABLE_ALL_WARNINGS)
  # GCC, Clang, and Intel seem to accept these
    list(APPEND REMORA_CXX_FLAGS "-Wall" "-Wextra" "-pedantic")
  if(NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    # ifort doesn't like -Wall
    list(APPEND REMORA_Fortran_FLAGS "-Wall")
  else()
    # Intel always reports some diagnostics we don't necessarily care about
    list(APPEND REMORA_CXX_FLAGS "-diag-disable:11074,11076")
  endif()
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 7.0)
    # Avoid notes about -faligned-new with GCC > 7
    list(APPEND REMORA_CXX_FLAGS "-faligned-new")
  endif()
endif()

# Add our extra flags according to language
separate_arguments(REMORA_CXX_FLAGS)
target_compile_options(${target} PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${REMORA_CXX_FLAGS}>)

separate_arguments(REMORA_Fortran_FLAGS)
target_compile_options(${target} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${REMORA_Fortran_FLAGS}>)

if(REMORA_ENABLE_CUDA)
  list(APPEND REMORA_CUDA_FLAGS "--expt-relaxed-constexpr")
  list(APPEND REMORA_CUDA_FLAGS "--expt-extended-lambda")
  list(APPEND REMORA_CUDA_FLAGS "--Wno-deprecated-gpu-targets")
  list(APPEND REMORA_CUDA_FLAGS "-m64")
  if(ENABLE_CUDA_FASTMATH)
    list(APPEND REMORA_CUDA_FLAGS "--use_fast_math")
  endif()
  separate_arguments(REMORA_CUDA_FLAGS)
  target_compile_options(${target} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:${REMORA_CUDA_FLAGS}>)
  # AMReX resolves the CUDA architectures itself (Tools/CMake/AMReXCUDAArchs.cmake) into
  # CMAKE_CUDA_ARCHITECTURES and exports the result as AMREX_CUDA_ARCHS. Build this target
  # for the same architectures as AMReX; if AMReX did not export them, the target inherits
  # CMAKE_CUDA_ARCHITECTURES.
  if(AMREX_CUDA_ARCHS)
    set_target_properties( ${target}
       PROPERTIES
       CUDA_ARCHITECTURES "${AMREX_CUDA_ARCHS}"
       )
  endif()
  set_target_properties(
    ${target} PROPERTIES
    LANGUAGE CUDA
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_RESOLVE_DEVICE_SYMBOLS ON)
endif()
endfunction()
