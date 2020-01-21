#!/bin/bash

np=4;
#as=8;
#res=80;
correction=0;

if [ $as -le 9 ]
then
	as=0$as
else
	as=$as
fi
echo $as

if [ $correction -eq 0 ]
then
	mpirun -n ${np} ../../../../parrot-opt -i 2advection.i resolution=${res} adapSteps=${as} UserObjects/active='soln'
fi
if [ $correction -eq 1 ]
then
	mpirun -n ${np} ../../../../parrot-opt -i 2advection.i resolution=${res} adapSteps=${as} Problem/operator_userobject='storeOperatorsUO'
fi
if [ $correction -eq 2 ]
then
	mpirun -n ${np} ../../../../parrot-opt -i 1diffusion.i resolution=${res} adapSteps=${as}
fi
