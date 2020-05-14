#!/bin/bash

# echo $refString

if [ $problemType -eq 0 ]
then
# mpirun -n ${np} ../../../parrot-opt -i 1diffusion.i mRefType=${refinementType} mRefLev=${refLev} MeshModifiers/my/refinements='14 0'
mpirun -n ${np} ../../../parrot-opt -i 1diffusion.i mRefType=0 mRefLev=0 MeshModifiers/my/refinements='14 0'
exit 
mpirun -n ${np} ../../../parrot-opt -i 1diffusion.i mRefType=0 mRefLev=1 MeshModifiers/my/refinements='15 0'
mpirun -n ${np} ../../../parrot-opt -i 1diffusion.i mRefType=0 mRefLev=2 MeshModifiers/my/refinements='16 0'
mpirun -n ${np} ../../../parrot-opt -i 1diffusion.i mRefType=0 mRefLev=3 MeshModifiers/my/refinements='17 0'
mpirun -n ${np} ../../../parrot-opt -i 1diffusion.i mRefType=0 mRefLev=4 MeshModifiers/my/refinements='18 0'
mpirun -n ${np} ../../../parrot-opt -i 1diffusion.i mRefType=1 mRefLev=0 MeshModifiers/my/refinements='1 1 1 1 1 1 1 1 1 1 1 1 1 1'
mpirun -n ${np} ../../../parrot-opt -i 1diffusion.i mRefType=1 mRefLev=1 MeshModifiers/my/refinements='1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0'
mpirun -n ${np} ../../../parrot-opt -i 1diffusion.i mRefType=1 mRefLev=2 MeshModifiers/my/refinements='1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1'
mpirun -n ${np} ../../../parrot-opt -i 1diffusion.i mRefType=1 mRefLev=3 MeshModifiers/my/refinements='1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0'
mpirun -n ${np} ../../../parrot-opt -i 1diffusion.i mRefType=1 mRefLev=4 MeshModifiers/my/refinements='1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1'
fi

if [ $problemType -eq 1 ]
then
	mpirun -n ${np} ../../../parrot-opt -i 2advectionM.i mRefType=${refinementType} mRefLev=${refLev} MeshModifiers/my/refinements=${refString}
fi

