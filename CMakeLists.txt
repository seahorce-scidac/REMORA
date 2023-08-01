############################ BASE ######################################

cmake_minimum_required (VERSION 3.14 FATAL_ERROR)
project(ROMSX CXX C)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake")
include(CMakePackageConfigHelpers)

########################## OPTIONS #####################################

#General options for all executables in the project
set(ROMSX_DIM "3" CACHE STRING "Number of physical dimensions")
option(ROMSX_ENABLE_DOCUMENTATION "Build documentation" OFF)
option(ROMSX_ENABLE_SALINITY      "Enable Salinity" ON)
option(ROMSX_ENABLE_ALL_WARNINGS  "Enable all compiler warnings" OFF)
option(ROMSX_ENABLE_TESTS         "Enable regression and unit tests" OFF)
option(ROMSX_ENABLE_NETCDF        "Enable NetCDF IO" OFF)
option(ROMSX_ENABLE_HDF5          "Enable HDF5 IO" ${ROMSX_ENABLE_NETCDF})
option(ROMSX_ENABLE_FCOMPARE "Enable building fcompare when not testing" OFF)
set(ROMSX_PRECISION "DOUBLE" CACHE STRING "Floating point precision SINGLE or DOUBLE")

#Options for performance
option(ROMSX_ENABLE_MPI    "Enable MPI"    OFF)
option(ROMSX_ENABLE_OPENMP "Enable OpenMP" OFF)
option(ROMSX_ENABLE_CUDA   "Enable CUDA"   OFF)
option(ROMSX_ENABLE_HIP    "Enable HIP"    OFF)
option(ROMSX_ENABLE_SYCL   "Enable SYCL"   OFF)

#Options for C++
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(ROMSX_ENABLE_CUDA)
  enable_language(CUDA)
  if(CMAKE_CUDA_COMPILER_VERSION VERSION_LESS "9.0")
    message(FATAL_ERROR "Your nvcc version is ${CMAKE_CUDA_COMPILER_VERSION} which is unsupported."
      "Please use CUDA toolkit version 9.0 or newer.")
  endif()
endif()

if(NOT ROMSX_DIM EQUAL 3)
  message(FATAL_ERROR "ROMSX is only supported in 3D.")
endif()

# Configure measuring code coverage in tests
option(CODECOVERAGE "Enable code coverage profiling" OFF)
if(CODECOVERAGE)
  # Only supports GNU
  if(NOT CMAKE_CXX_COMPILER_ID MATCHES GNU)
    message(WARNING "CODECOVERAGE is only support with GNU Compilers. The current C++ compiler is ${CMAKE_CXX_COMPILER_ID}")
  endif()
  if(NOT CMAKE_C_COMPILER_ID MATCHES GNU)
    message(WARNING "CODECOVERAGE is only support with GNU Compilers. The current C compiler is ${CMAKE_C_COMPILER_ID}")
  endif()

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")
endif()

########################### AMReX #####################################

set(AMREX_SUBMOD_LOCATION "${CMAKE_SOURCE_DIR}/Submodules/AMReX")
include(${CMAKE_SOURCE_DIR}/CMake/SetAmrexOptions.cmake)
list(APPEND CMAKE_MODULE_PATH "${AMREX_SUBMOD_LOCATION}/Tools/CMake")

########################### AMReX #####################################

add_subdirectory(${AMREX_SUBMOD_LOCATION})
include(${CMAKE_SOURCE_DIR}/CMake/SetAmrexCompileFlags.cmake)

########################## NETCDF ##################################

if(ROMSX_ENABLE_NETCDF)
  set(CMAKE_PREFIX_PATH ${NETCDF_DIR} ${CMAKE_PREFIX_PATH})
##  set(NETCDF_CXX "YES")
  find_package (NetCDF REQUIRED)
  if(NETCDF_FOUND)
     message(STATUS "Found NetCDF, NETCDF_DIR = ${NETCDF_DIR}")
  endif()
endif()

########################### ROMSX #####################################

if(ROMSX_ENABLE_MPI)
  find_package(MPI REQUIRED)
endif()

# General information about machine, compiler, and build type
message(STATUS "ROMSX Information:")
message(STATUS "CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
message(STATUS "CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION = ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

# Turn on rpath stuff
include(${CMAKE_SOURCE_DIR}/CMake/SetRpath.cmake)

#Build ROMSX executables and intermediate object library and link to amrex library
add_subdirectory(Exec)

if(ROMSX_ENABLE_TESTS)
  include(CTest)
  add_subdirectory(Tests)
endif()

if(ROMSX_ENABLE_DOCUMENTATION)
   add_subdirectory(Docs)
endif()
