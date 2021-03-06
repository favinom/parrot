#!/bin/bash

as=5;
us=1;
np=1;
mpirun -n ${np} ../../../parrot-opt -i 0refineBlock.i adapSteps=${as} unifSteps=${us}

for (( c=0; c<=as+1; c++ ))
do
rm refinedMesh_${us}_000$c.xdr
done

for (( c=1; c<=as+1; c++ ))
do
temp=`expr $c - 1`
mv refinedMesh_${us}_000${c}_mesh.xdr refinedMesh_${us}_000${temp}_mesh.xdr
done

mpirun -n ${np} ../../../parrot-opt -i 1diffusion.i adapSteps=${as} unifSteps=${us}
#mpirun -n ${np} ../../parrot-opt -i 2advection.i adapSteps=${as} unifSteps=${us}
#mpirun -n ${np} ../../parrot-opt -i 2advectionAFC.i adapSteps=${as} unifSteps=${us}

