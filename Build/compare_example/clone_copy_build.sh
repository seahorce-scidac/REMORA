#!/bin/bash

set -e

umask 002

git clone --recursive git@github.com:seahorce-scidac/ROMSX
git clone --recursive git@github.com:seahorce-scidac/COAWST

cd ROMSX/Build/compare_example

cp upwelling.h ../../../COAWST/ROMS/Include/upwelling.h

cp roms_upwelling.in ../../../COAWST/ROMS/External/roms_upwelling.in

cp ana_grid.h ../../../COAWST/ROMS/Functionals/ana_grid.h
cp ana_initial.h ../../../COAWST/ROMS/Functionals/ana_initial.h
cp *.F ../../../COAWST/ROMS/Nonlinear/
cp *.mk ../../../COAWST/Compilers
cp coawst.bash ../../../COAWST

cd ../
if [ "$NERSC_HOST" == "cori" ]
then
source cori-env.sh
elif [ "$NERSC_HOST" == "perlmutter" ]
then
source saul-env.sh
else
export NETCDF_INCDIR=${NETCDF_DIR}/include
export NETCDF_LIBDIR=${NETCDF_DIR}/lib
fi

echo "We will use these paths for netcdf:"
echo "${NETCDF_DIR}"
echo "${NETCDF_LIBDIR}"
echo "${NETCDF_INCDIR}"

cd ../Exec/Upwelling
nice make -j16 USE_NETCDF=TRUE DEBUG=TRUE
cd ../../../

cd COAWST

#activeList = {
#  ["Nsight-Compute"] = "2022.1.1",
#  ["Nsight-Systems"] = "2022.2.1",
#  ["PrgEnv-gnu"] = "8.3.3",
#  ["cpe"] = "22.11",
#  ["cray-dsmml"] = "0.2.2",
#  ["cray-hdf5"] = "1.12.2.1",
#  ["cray-libsci"] = "22.11.1.2",
#  ["cray-mpich"] = "8.1.22",
#  ["cray-netcdf"] = "4.9.0.1",
#  ["craype"] = "2.7.19",
#  ["craype-accel-nvidia80"] = "",
#  ["craype-network-ofi"] = "",
#  ["craype-x86-milan"] = "",
#  ["cudatoolkit"] = "11.7",
#  ["gcc"] = "11.2.0",
#  ["gpu"] = "1.0",
#  ["libfabric"] = "1.15.2.0",
#  ["perftools-base"] = "22.09.0",
#  ["xalt"] = "2.10.2",
#  ["xpmem"] = "2.5.2-2.4_3.20__gd0f7936.shasta",
#}

./coawst.bash -j 4
