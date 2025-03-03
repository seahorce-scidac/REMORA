############################ BASE ######################################

cmake_minimum_required (VERSION 3.14 FATAL_ERROR)
project(REMORA CXX C)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake")
include(CMakePackageConfigHelpers)

########################## OPTIONS #####################################

#General options for all executables in the project
set(REMORA_DIM "3" CACHE STRING "Number of physical dimensions")
option(REMORA_ENABLE_DOCUMENTATION "Build documentation" OFF)
option(REMORA_ENABLE_ALL_WARNINGS  "Enable all compiler warnings" OFF)
option(REMORA_ENABLE_TESTS         "Enable regression and unit tests" OFF)
option(REMORA_ENABLE_PNETCDF        "Enable PNetCDF IO" OFF)
option(REMORA_ENABLE_FCOMPARE "Enable building fcompare when not testing" OFF)
set(REMORA_PRECISION "DOUBLE" CACHE STRING "Floating point precision SINGLE or DOUBLE")

#Options for performance
option(REMORA_ENABLE_MPI    "Enable MPI"    OFF)
option(REMORA_ENABLE_OPENMP "Enable OpenMP" OFF)
option(REMORA_ENABLE_CUDA   "Enable CUDA"   OFF)
option(REMORA_ENABLE_HIP    "Enable HIP"    OFF)
option(REMORA_ENABLE_SYCL   "Enable SYCL"   OFF)
option(REMORA_ENABLE_MOAB   "Enable MOAB"   OFF)

#Options for C++
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(REMORA_ENABLE_CUDA)
  enable_language(CUDA)
  if(CMAKE_CUDA_COMPILER_VERSION VERSION_LESS "9.0")
    message(FATAL_ERROR "Your nvcc version is ${CMAKE_CUDA_COMPILER_VERSION} which is unsupported."
      "Please use CUDA toolkit version 9.0 or newer.")
  endif()
endif()

if(NOT REMORA_DIM EQUAL 3)
  message(FATAL_ERROR "REMORA is only supported in 3D.")
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

########################## PNETCDF ##################################

if(REMORA_ENABLE_PNETCDF)
  find_package (PNetCDF REQUIRED)
  if(PNETCDF_FOUND)
     message(STATUS "Found PNetCDF, PNETCDF_DIR = ${PNETCDF_DIR}")
     message(STATUS " PNETCDF_INCLUDE_DIRS = ${PNETCDF_INCLUDE_DIRS}")
     message(STATUS " PNETCDF_VERSION = ${PNETCDF_VERSION}")
     message(STATUS " PNETCDF_LIBRARY_DIRS = ${PNETCDF_LIBRARY_DIRS}")
     message(STATUS " PNETCDF_ROOT_DIR = ${PNETCDF_ROOT_DIR}")
  endif()
endif()

if ( REMORA_ENABLE_MOAB )
  find_package( MOAB REQUIRED )
  # set(REMORA_USE_MOAB ON)
  set( REMORA_DEFINES " -DREMORA_USE_MOAB ${REMORA_DEFINES}" )
  message(STATUS "MOAB cmake found at  ${MOAB_DIR} ")
  message(STATUS "MOAB_LIBRARIES:   ${MOAB_LIBRARIES} ")
  message(STATUS "MOAB_INCLUDE_DIRS:   ${MOAB_INCLUDE_DIRS} ")
endif( REMORA_ENABLE_MOAB )

########################### REMORA #####################################

if(REMORA_ENABLE_MPI)
  find_package(MPI REQUIRED)
endif()

# General information about machine, compiler, and build type
message(STATUS "REMORA Information:")
message(STATUS "CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
message(STATUS "CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION = ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

# Turn on rpath stuff
include(${CMAKE_SOURCE_DIR}/CMake/SetRpath.cmake)

#Build REMORA executables and intermediate object library and link to amrex library
add_subdirectory(Exec)

if(REMORA_ENABLE_TESTS)
  include(CTest)
  add_subdirectory(Tests)
endif()

if(REMORA_ENABLE_DOCUMENTATION)
   add_subdirectory(Docs)
endif()
