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

mpirun -n ${np} ../../../../parrot-opt -i 3advectionM.i resolution=${res} unifSteps=${us} adaptSteps=${as} Problem/operator_userobject='storeOperatorsUO' UserObjects/active='soln assembleVolumeVectors MassAssembly storeOperatorsUO' # assembleVolumeVectors2'
