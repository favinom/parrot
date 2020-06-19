#!/bin/bash

if [ $problemType -eq 0 ]
then
$runString scripts/0generate.i mRes=${res} mResName=${resName}
#$runString scripts/2diffusion.i mRes=${res} mResName=${resName} mRefLev=${amr} mRefLevName=${amrName} mUmr=${umr} # MeshModifiers/my/refinements=\'${val}\'
$runString scripts/1refine.i mRes=${res} mResName=${resName} mRefLev=${amr}
$runString scripts/testF.i mRes=${res} mResName=${resName} mRefLev=${amr} mRefLevName=${amrName} mUmr=${umr} # MeshModifiers/my/refinements=\'${val}\'
fi

if [ $problemType -eq 1 ]
then
#$runString scripts/0generate.i mRes=${res} mResName=${resName}
#$runString scripts/1refine.i mRes=${res} mResName=${resName} mRefLev=${amr}
$runString scripts/2advectionM.i mRes=${res} mRefLev=${amr} mRefLevName=${amrName} mUmr=${umr} phiIn=1.0 mResName=${resName}
fi

