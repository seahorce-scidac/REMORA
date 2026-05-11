############################ BASE ######################################

cmake_minimum_required (VERSION 3.14 FATAL_ERROR)
if(REMORA_ENABLE_FUNWAVE_FORT)
  project(REMORA CXX C Fortran)
else()
  project(REMORA CXX C)
endif()
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMake")
include(CMakePackageConfigHelpers)

########################## OPTIONS #####################################

#General options for all executables in the project
set(REMORA_DIM "3" CACHE STRING "Number of physical dimensions")
option(REMORA_ENABLE_DOCUMENTATION "Build documentation" OFF)
option(REMORA_ENABLE_ALL_WARNINGS  "Enable all compiler warnings" OFF)
option(REMORA_ENABLE_TESTS         "Enable regression and unit tests" OFF)
option(REMORA_ENABLE_PNETCDF        "Enable PNetCDF IO" OFF)
option(REMORA_USE_INTERNAL_AMREX   "Add AMReX as subproject" ON)

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(REMORA_BUILD_EXECUTABLES_DEFAULT ON)
else()
  set(REMORA_BUILD_EXECUTABLES_DEFAULT OFF)
endif()
option(REMORA_BUILD_EXECUTABLES "Build REMORA example executables" ${REMORA_BUILD_EXECUTABLES_DEFAULT})
set(REMORA_LIBRARY_PROBLEM "BlankProblem" CACHE STRING
    "Exec problem directory used to provide prob.cpp when REMORA_BUILD_EXECUTABLES=OFF")
set_property(CACHE REMORA_LIBRARY_PROBLEM PROPERTY STRINGS
    BlankProblem Advection Seamount DoublyPeriodic Upwelling Channel_Test
    DoubleGyre IdealMiniGrid BoundaryLayer DogboneAnalytic)

if(DEFINED REMORA_ENABLE_NETCDF)
  message(FATAL_ERROR "REMORA_ENABLE_NETCDF is not supported. Only REMORA_ENABLE_PNETCDF is supported.")
endif()
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

if (REMORA_USE_INTERNAL_AMREX)
  set(AMREX_SUBMOD_LOCATION "${PROJECT_SOURCE_DIR}/Submodules/AMReX")
  include(${PROJECT_SOURCE_DIR}/CMake/SetAmrexOptions.cmake)
  list(APPEND CMAKE_MODULE_PATH "${AMREX_SUBMOD_LOCATION}/Tools/CMake")
  add_subdirectory(${AMREX_SUBMOD_LOCATION})
  set(REMORA_AMREX_TARGET amrex)
else()
  if (TARGET AMReX::amrex OR TARGET amrex)
    if (TARGET AMReX::amrex)
      set(REMORA_AMREX_TARGET AMReX::amrex)
    else()
      set(REMORA_AMREX_TARGET amrex)
    endif()
  else()
    set(CMAKE_PREFIX_PATH ${AMReX_DIR} ${CMAKE_PREFIX_PATH})
    list(APPEND AMREX_COMPONENTS "3D" "DOUBLE")
    if (REMORA_ENABLE_MPI)
      list(APPEND AMREX_COMPONENTS "MPI")
    endif()
    if (REMORA_ENABLE_OPENMP)
      list(APPEND AMREX_COMPONENTS "OMP")
    endif()
    if (REMORA_ENABLE_CUDA)
      list(APPEND AMREX_COMPONENTS "CUDA")
    endif()
    if (REMORA_ENABLE_SYCL)
      list(APPEND AMREX_COMPONENTS "SYCL")
    endif()
    if (REMORA_ENABLE_HIP)
      list(APPEND AMREX_COMPONENTS "HIP")
    endif()
    separate_arguments(AMREX_COMPONENTS)
    find_package(AMReX CONFIG REQUIRED COMPONENTS ${AMREX_COMPONENTS})
    set(REMORA_AMREX_TARGET AMReX::amrex)
  endif()
endif()

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
include(${PROJECT_SOURCE_DIR}/CMake/SetRpath.cmake)

#Build REMORA executables and intermediate object library and link to amrex library
add_subdirectory(Exec)

if(REMORA_ENABLE_TESTS)
  include(CTest)
  add_subdirectory(Tests)
endif()

if(REMORA_ENABLE_DOCUMENTATION)
   add_subdirectory(Docs)
endif()

# Installation rules
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

# Create non-object library for use as external target
add_library(remora_api)
if (BUILD_SHARED_LIBS)
  set_target_properties(remora_api PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()
if (REMORA_ENABLE_CUDA)
   setup_target_for_cuda_compilation(remora_api)
endif ()
target_link_libraries(remora_api PUBLIC remora_srclib)
add_library(${PROJECT_NAME}::remora_api ALIAS remora_srclib)

# Collect all headers and make them installable with the target
file(GLOB_RECURSE REMORA_INCLUDES "${PROJECT_SOURCE_DIR}/Source/*.H")

set_target_properties(
  remora_srclib PROPERTIES PUBLIC_HEADER "${REMORA_INCLUDES}")

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  # Install REMORA
  install(
    TARGETS remora_api remora_srclib
    EXPORT ${PROJECT_NAME}Targets
    RUNTIME       DESTINATION bin
    ARCHIVE       DESTINATION lib
    LIBRARY       DESTINATION lib
    INCLUDES      DESTINATION include
    PUBLIC_HEADER DESTINATION include
    )

  # Install all headers in include directories
  #install(
  #  DIRECTORY ${REMORA_INCLUDE_DIRECTORIES}
  #  DESTINATION include
  #  FILES_MATCHING PATTERN "*.H")

  # Make REMORA discoverable using `find_package`
  install(
    EXPORT ${PROJECT_NAME}Targets
    NAMESPACE ${PROJECT_NAME}::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
    )
  configure_package_config_file(
    CMake/${PROJECT_NAME}Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
    )
  install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
    )
endif()
