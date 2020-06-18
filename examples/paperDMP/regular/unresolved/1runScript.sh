#!/bin/bash

if [ $problemType -eq 0 ]
then
$runString scripts/1diffusion.i mRes=${res} mResName=${resName} mRefLev=${amr} mRefLevName=${amrName} mUmr=${umr} # MeshModifiers/my/refinements=\'${val}\'
fi

if [ $problemType -eq 1 ]
then
$runString scripts/2advection.i mRes=${res} mResName=${resName} mRefLev=${amr} mRefLevName=${amrName} mUmr=${umr} # MeshModifiers/my/refinements=\'${val}\'
fi
