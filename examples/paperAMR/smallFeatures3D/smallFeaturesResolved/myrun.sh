#!/bin/bash

np=4;
createmesh=1;
correction=0;


#mpirun -n ${np} ../../../parrot-opt -i 1diffusion.i
#exit 1

if [ $correction -eq 0 ]
then
	mpirun -n ${np} ../../../parrot-opt -i 2advection.i adapSteps=${as}
fi
if [ $correction -eq 1 ]
then
	mpirun -n ${np} ../../../parrot-opt -i 2advectionCorrection.i adapSteps=${as}
fi
