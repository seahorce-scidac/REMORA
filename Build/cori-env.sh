module unload PrgEnv-intel
module load PrgEnv-gnu
module load cray-hdf5-parallel
module load cray-netcdf-hdf5parallel
module load ncview
module load openmpi

export NETCDF_INCDIR=${NETCDF_DIR}/include
export NETCDF_LIBDIR=${NETCDF_DIR}/lib
