#!/bin/bash

# number of processors
np=2;
#background elements, THIS NUMBER HAS TO BE DIVISIBLE BY 8
be=600;
# element per fracture
fe=4;
# starting level
res=`expr ${be} + 3 \\* ${fe}`
#echo $res

# number of adaptive steps
as=1;
us=0;

rm -rf uniform*

mpirun -n ${np} ../../parrot-opt -i 0generateMesh.i resolution=${res}

rm uniform_${res}_0000.xda
rm uniform_${res}_0001.xda
rm uniform_${res}_0000_mesh.xda
mv uniform_${res}_0001_mesh.xda uniform_${res}_mesh.xda

./separateFile.sh uniform_${res}_mesh.xda ${res}

commandString='faimesh('${be}','${fe}')'

#/Users/mariagiuseppinanestola/Desktop/MATLAB_R2014b.app/bin/matlab -nodesktop -nosplash -r "faimesh; quit"
/Applications/MATLAB_R2018b.app/bin/matlab -nodesktop -nosplash -r ${commandString}

rm uniform_${res}_mesh.xda

./uniteFile.sh adapted_${res}_mesh.xda

mpirun -n ${np} ../../parrot-opt -i 1refineBlock.i resolution=${res} unifSteps=${us} adaptSteps=${as}

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

#mpirun -n ${np} ../../parrot-opt -i 2diffusion.i resolution=${res} unifSteps=${us} adaptSteps=${as}
#mpirun -n ${np} ../../parrot-opt -i 2advection.i resolution=${res} unifSteps=${us} adaptSteps=${as}
mpirun -n ${np} ../../parrot-opt -i 2advection_afc.i resolution=${res} unifSteps=${us} adaptSteps=${as}