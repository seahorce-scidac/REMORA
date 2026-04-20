#!/bin/bash
CXX_cmake=mpicxx
C_cmake=mpicc
FC_cmake=mpifort
if [ "$NERSC_HOST" == "perlmutter" ]
then
# See https://docs.nersc.gov/development/compilers/wrappers/#hpe-cray-compiler-wrappers
     source ../Build/saul-env.sh
     CXX_cmake=CC
     C_cmake=cc
     FC_cmake=ftn
fi
# Example CMake config script for an OSX laptop with OpenMPI

cmake -DCMAKE_INSTALL_PREFIX:PATH=./install \
      -DCMAKE_CXX_COMPILER:STRING=${CXX_cmake} \
      -DCMAKE_C_COMPILER:STRING=${C_cmake} \
      -DCMAKE_Fortran_COMPILER:STRING=${FC_cmake} \
      -DMPIEXEC_PREFLAGS:STRING=--oversubscribe \
      -DCMAKE_BUILD_TYPE:STRING=Release \
      -DREMORA_DIM:STRING=3 \
      -DREMORA_ENABLE_MPI:BOOL=ON \
      -DREMORA_ENABLE_TESTS:BOOL=ON \
      -DREMORA_ENABLE_FCOMPARE:BOOL=ON \
      -DREMORA_ENABLE_DOCUMENTATION:BOOL=OFF \
      -DREMORA_ENABLE_PNETCDF:BOOL=ON \
      -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON \
      .. && make -j8 && make install
