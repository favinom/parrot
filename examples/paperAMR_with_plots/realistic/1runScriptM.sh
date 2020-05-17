#!/bin/bash

if [ $problemType -eq 0 ]
then
$runline scripts/1diffusion.i mRefLev=${amr} mRefLevName=${amrName} mUmr=${umr} # MeshModifiers/my/refinements=\'${val}\'
fi

if [ $problemType -eq 1 ]
then
$runline scripts/2advectionM.i mRefLev=${amr} mRefLevName=${amrName} mUmr=${umr} phiIn=0.00001157407
fi

