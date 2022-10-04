#!/bin/bash
export inFilename=roms_upwelling.in
export inputFilename=inputs
export clearInputs=FALSE
if [ ${clearInputs} != FALSE ]; then
   /dev/null > inputs;
fi

export varNameList="S0 T0"
for varName in ${varNameList}; do
    grep ${varName} ${inFilename} | grep == | sed s/==/=/g | sed s/"!"/"#"/g | sed s/${varName}/prob.${varName}/g >> inputs
done
