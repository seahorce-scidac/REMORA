module load PrgEnv-gnu
module load cray-parallel-netcdf
# module load cudatoolkit   # GPU builds
export PKG_CONFIG_PATH=${MPICH_DIR}/lib/pkgconfig:$PKG_CONFIG_PATH
