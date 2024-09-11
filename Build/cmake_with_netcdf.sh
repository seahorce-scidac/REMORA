CXX_cmake=mpicxx
C_cmake=mpicx
FC_cmake=mpifort
if [ "$NERSC_HOST" == "perlmutter" ]
then
# See https://docs.nersc.gov/development/compilers/wrappers/#hpe-cray-compiler-wrappers
     source saul-env.sh
     CXX_cmake=CC
     C_cmake=cc
     FC_cmake=ftn
     #     cp /opt/cray/pe/netcdf-hdf5parallel/4.9.0.9/gnu/12.3/lib/pkgconfig/netcdf-cxx4_parallel.pc netcdf.pc
     cp ${CRAY_NETCDF_HDF5PARALLEL_PREFIX}/lib/pkgconfig/${PE_NETCDF_HDF5PARALLEL_CXX_PKGCONFIG_LIBS}.pc netcdf.pc
     sed -i s/netcdf-cxx4/netcdf/g netcdf.pc
     export PKG_CONFIG_PATH=$(pwd):${PE_GNU_FIXED_PKGCONFIG_PATH}:${PKG_CONFIG_PATH}
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
      -DREMORA_ENABLE_NETCDF:BOOL=ON \
      -DREMORA_ENABLE_HDF5:BOOL=ON \
      -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON \
      .. && make -j8
