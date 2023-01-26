#!/bin/bash

set -e

umask 002

git clone --recursive git@github.com:seahorce-scidac/ROMSX
git clone --recursive git@github.com:seahorce-scidac/COAWST

cd ROMSX/Build/compare_example
diff upwelling.h ../../../COAWST/ROMS/Include/upwelling.h
cp upwelling.h ../../../COAWST/ROMS/Include/upwelling.h
diff roms_upwelling.in ../../../COAWST/ROMS/External/roms_upwelling.in
cp roms_upwelling.in ../../../COAWST/ROMS/External/roms_upwelling.in
diff ana_grid.h ../../../COAWST/ROMS/Functionals/ana_grid.h
cp ana_grid.h ../../../COAWST/ROMS/Functionals/ana_grid.h
cp *.F ../../../COAWST/ROMS/Nonlinear/
cp *.mk ../../../COAWST/Compilers

cd ../
source saul-env.sh
cd ../Exec/Upwelling
nice make -j16 USE_NETCDF=TRUE USE_CUDA=TRUE
cd ../../../

cd COAWST
sed -i s/'\/global\/homes\/h\/hetland\/COAWST'/'$(pwd)'/g coawst.bash
git checkout ROMS/Include/upwelling.h

module swap cray-hdf5-parallel cray-hdf5
module swap cray-netcdf-hdf5parallel cray-netcdf

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

./coawst.bash -j 20
