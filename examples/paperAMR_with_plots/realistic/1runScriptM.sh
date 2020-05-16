#!/bin/bash

if [ $problemType -eq 0 ]
then
$runline scripts/1diffusion.i mRefLev=${refLev} mRefLevName=${refLevName} # MeshModifiers/my/refinements=\'${val}\'
fi

if [ $problemType -eq 1 ]
then
$runline scripts/2advectionM.i mRefLev=${refLev} mRefLevName=${refLevName}
fi

