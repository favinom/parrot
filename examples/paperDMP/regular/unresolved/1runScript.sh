#!/bin/bash

if [ $problemType -eq 0 ]
then
$runString scripts/0generate.i mRes=${res} mResName=${resName}
$runString scripts/2diffusion.i mRes=${res} mResName=${resName} mRefLev=${amr} mRefLevName=${amrName} mUmr=${umr} # MeshModifiers/my/refinements=\'${val}\'
fi

if [ $problemType -eq 1 ]
then
$runString scripts/2advectionM.i mRefLev=${amr} mRefLevName=${amrName} mUmr=${umr} phiIn=1.0
fi

