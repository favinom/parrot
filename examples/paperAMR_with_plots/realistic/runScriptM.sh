#!/bin/bash

np=4;

if [ $as -le 9 ]
then
	as=0$as
else
	as=$as
fi

mpirun -n ${np} ../../../parrot-opt -i 1diffusion.i resolution=${res} adapSteps=${as}
#  Problem/operator_userobject='storeOperatorsUO' UserObjects/active='soln assembleVolumeVectors MassAssembly storeOperatorsUO'

# mpirun -n ${np} ../../../../parrot-opt -i 2advectionM.i resolution=${res} adapSteps=${as} Problem/operator_userobject='storeOperatorsUO' UserObjects/active='soln assembleVolumeVectors MassAssembly storeOperatorsUO'
