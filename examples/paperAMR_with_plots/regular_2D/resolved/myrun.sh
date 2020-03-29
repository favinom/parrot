#!/bin/bash

# number of processors
np=4;
#background elements, THIS NUMBER HAS TO BE DIVISIBLE BY 8
#be=200;
# element per fracture
#fe=8;
# starting level
res=`expr ${be} + 3 \\* ${fe}`
#echo $res

# number of adaptive steps
as=0;
us=0;

createmesh=1;
correction=4;

if [ $createmesh -eq 1 ]
then
	rm -rf uniform*

	mpirun -n ${np} ../../../../parrot-opt -i 0generateMesh.i resolution=${res}

	rm uniform_${res}_0000.xda
	rm uniform_${res}_0001.xda
	rm uniform_${res}_0000_mesh.xda
	mv uniform_${res}_0001_mesh.xda uniform_${res}_mesh.xda

	./separateFile.sh uniform_${res}_mesh.xda ${res}

	commandString='faimesh('${be}','${fe}')'

   /Users/mariagiuseppinanestola/Desktop/MATLAB_R2014b.app/bin/matlab -nodesktop -nosplash -r ${commandString}
   #/Applications/MATLAB_R2018b.app/bin/matlab -nodesktop -nosplash -r ${commandString}

	rm uniform_${res}_mesh.xda

	./uniteFile.sh adapted_${res}_mesh.xda

	mpirun -n ${np} ../../../../parrot-opt -i 1refineBlock.i resolution=${res} unifSteps=${us} adaptSteps=${as}

	for (( c=0; c<=as+1; c++ ))
	do
		rm refined_${res}_${us}_000$c.xdr
	done

	for (( c=1; c<=as+1; c++ ))
	do
		temp=`expr $c - 1`
		mv refined_${res}_${us}_000${c}_mesh.xdr refined_${res}_${us}_000${temp}_mesh.xdr
	done

	for (( c=0; c<as; c++ ))
	do
		rm refined_${res}_${us}_000${c}_mesh.xdr
	done
fi

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
