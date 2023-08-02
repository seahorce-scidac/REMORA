#!/bin/bash

# Example CMake config script for an OSX laptop with OpenMPI

cmake -DCMAKE_INSTALL_PREFIX:PATH=./install \
      -DCMAKE_CXX_COMPILER:STRING=mpicxx \
      -DCMAKE_C_COMPILER:STRING=mpicc \
      -DCMAKE_Fortran_COMPILER:STRING=mpifort \
      -DMPIEXEC_PREFLAGS:STRING=--oversubscribe \
      -DCMAKE_BUILD_TYPE:STRING=Release \
      -DROMSX_DIM:STRING=3 \
      -DROMSX_ENABLE_MPI:BOOL=ON \
      -DROMSX_ENABLE_TESTS:BOOL=ON \
      -DROMSX_ENABLE_FCOMPARE:BOOL=ON \
      -DROMSX_ENABLE_DOCUMENTATION:BOOL=OFF \
      -DROMSX_ENABLE_SALINITY:BOOL=ON \
      -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON \
      .. && make -j8
