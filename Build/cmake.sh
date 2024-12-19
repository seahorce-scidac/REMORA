#!/bin/bash

# Example CMake config script for an OSX laptop with OpenMPI

cmake -DCMAKE_INSTALL_PREFIX:PATH=./install \
      -DCMAKE_CXX_COMPILER:STRING=mpicxx \
      -DMPIEXEC_PREFLAGS:STRING=--oversubscribe \
      -DCMAKE_BUILD_TYPE:STRING=Release \
      -DREMORA_DIM:STRING=3 \
      -DREMORA_ENABLE_MPI:BOOL=ON \
      -DREMORA_ENABLE_TESTS:BOOL=ON \
      -DREMORA_ENABLE_FCOMPARE:BOOL=ON \
      -DREMORA_ENABLE_DOCUMENTATION:BOOL=OFF \
      -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON \
      .. && make -j8
