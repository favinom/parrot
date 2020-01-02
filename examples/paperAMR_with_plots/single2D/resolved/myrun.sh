#!/bin/bash

as=0;
us=2;
or=5;
np=6;
createmesh=1;
correction=1;

if [ $createmesh -eq 1 ]
then
   ./myclean.sh
	mpirun -n ${np} ../../../../parrot-opt -i 0refineBlock.i typem=${or} adapSteps=${as} unifSteps=${us}
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
	mpirun -n ${np} ../../../../parrot-opt -i 2advection.i typem=${or} adapSteps=${as} unifSteps=${us}
fi
if [ $correction -eq 1 ]
then
    mpirun -n ${np} ../../../../parrot-opt -i 2advectionCorrection.i typem=${or} adapSteps=${as} unifSteps=${us}
    mpirun -n ${np} ../../../../parrot-opt -i plotLine_master_P.i typem=${or} adapSteps=${as} unifSteps=${us}
    mpirun -n ${np} ../../../../parrot-opt -i plotLine_master_CM.i typem=${or} adapSteps=${as} unifSteps=${us}
    mpirun -n ${np} ../../../../parrot-opt -i plotLine_master_CM2.i typem=${or} adapSteps=${as} unifSteps=${us}
fi

