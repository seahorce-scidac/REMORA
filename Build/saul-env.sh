module load PrgEnv-gnu
module load cray-hdf5-parallel
module load cray-netcdf-hdf5parallel
#module load ncview
#module load openmpi

export NETCDF_INCDIR=${NETCDF_DIR}/include
export NETCDF_LIBDIR=${NETCDF_DIR}/lib

module load cray-mpich
export PATH=${MPICH_DIR}/bin:$PATH
export LD_LIBRARY_PATH=${MPICH_DIR}/lib:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=/opt/cray/pe/mpich/8.1.17/ofi/gnu/9.1/lib/pkgconfig/:$PKG_CONFIG_PATH
