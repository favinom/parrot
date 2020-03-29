#!/bin/bash
as=3;
us=1;
np=4;
# this is type mesh: 0 old style, 1 the better one
tm=1;
or=07;

createmesh=1;
correction=0;

if [ $createmesh -eq 1 ]
then
	./myclean.sh
	mpirun -n ${np} ../../../../parrot-opt -i 0refineBlock.i adapSteps=${as} unifSteps=${us} tipeMesh=${tm} origRef=${or}
	for (( c=0; c<=as+1; c++ ))
	do
		rm refinedMesh_${us}_000$c.xdr
	done

	for (( c=1; c<=as+1; c++ ))
	do
		temp=`expr $c - 1`
		mv refinedMesh_${us}_000${c}_mesh.xdr refinedMesh_${us}_000${temp}_mesh.xdr
	done
fi

echo $as
#mpirun -n ${np} ../../../../parrot-opt -i 1diffusion.i adapSteps=${as} unifSteps=${us}
if [ $correction -eq 0 ]
then
	mpirun -n ${np} ../../../../parrot-opt -i 2advection.i adapSteps=${as} unifSteps=${us}
fi
if [ $correction -eq 1 ]
then
	mpirun -n ${np} ../../../../parrot-opt -i 2advectionCorrection.i adapSteps=${as} unifSteps=${us}
fi

