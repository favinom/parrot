#!/bin/bash

# number of processors
np=4;
#resolution
res=`expr ${be} + 3 \\* ${fe}`
#echo $res

# number of adaptive steps
as=0;
us=0;

createmesh=1;
correction=0;

if [ $correction -eq 0 ]
then
	#mpirun -n ${np} ../../../parrot-opt -i 2advection.i resolution=${res} unifSteps=${us} adaptSteps=${as}
	mpirun -n ${np} ../../../../parrot-opt -i 3advection.i resolution=${res} unifSteps=${us} adaptSteps=${as} UserObjects/active='soln'
fi
if [ $correction -eq 1 ]
then
	mpirun -n ${np} ../../../../parrot-opt -i 3advection.i resolution=${res} unifSteps=${us} adaptSteps=${as} Problem/operator_userobject='storeOperatorsUO'
fi
if [ $correction -eq 2 ]
then
	mpirun -n ${np} ../../../parrot-opt -i 2diffusion.i resolution=${res} unifSteps=${us} adaptSteps=${as}
fi
