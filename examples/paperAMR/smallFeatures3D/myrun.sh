#!/bin/bash

as=1;
np=4;
createmesh=1;
correction=1;

if [ $createmesh -eq 1 ]
then
	./myclean.sh
	mpirun -n ${np} ../../../parrot-opt -i 0refineBlock.i adapSteps=${as}
	for (( c=0; c<=as+1; c++ ))
	do
		if [ $c -le 9 ]
		then
			v=0$c
		else
			v=$c
		fi
		rm refinedMesh_00$v.xdr
	done

	for (( c=1; c<=as+1; c++ ))
	do
		temp=`expr $c - 1`
		if [ $c -le 9 ]
		then
			v1=0$c
		else
			v1=$c
		fi
		if [ $temp -le 9 ]
		then
			v2=0$temp
		else
			v2=$temp
		fi
		mv refinedMesh_00${v1}_mesh.xdr refinedMesh_00${v2}_mesh.xdr
	done
fi

if [ $as -le 9 ]
then
	as=0$as
else
	as=$as
fi
echo $as

#mpirun -n ${np} ../../../parrot-opt -i 1diffusion.i adapSteps=${as}
#exit 1

if [ $correction -eq 0 ]
then
	mpirun -n ${np} ../../../parrot-opt -i 2advection.i adapSteps=${as}
fi
if [ $correction -eq 1 ]
then
	mpirun -n ${np} ../../../parrot-opt -i 2advectionCorrection.i adapSteps=${as}
fi
