#!/bin/bash -l
#BATCH --job-name=contact
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nestom@usi.ch
#SBATCH --time=24:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=10
#SBATCH  --partition=fat
#SBATCH --mem=100000

typem=0;

ol=1; #2

np=20;

as=0;

us=3;

correction=0;

createmesh=1;

postprocessor=1;


if [ $createmesh -eq 1 ]
then
    mpirun -n ${np} ../../../../parrot-opt -i 0refineBlock.i typeMesh=${typem} origLevel=${ol} adapSteps=${as} Uref=${us}
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

if [ $correction -eq 0 ]
then
	mpirun -n ${np} ../../../../parrot-opt -i 2advection.i typeMesh=${typem} origLevel=${ol} adapSteps=${as} Uref=${us}
fi
if [ $correction -eq 1 ]
then
	srun -n ${np} ../../../../parrot-opt -i 2advectionCorrection.i typeMesh=${typem} origLevel=${ol} adapSteps=${as} Uref=${us}
fi
if [ $postprocessor -eq 1 ]
then
    mpirun -n ${np} ../../../../parrot-opt -i plotLine_master_P.i typeMesh=${typem} origLevel=${ol} adapSteps=${as} Uref=${us}
    mpirun -n ${np} ../../../../parrot-opt -i plotLine_master_P_2.i typeMesh=${typem} origLevel=${ol} adapSteps=${as} Uref=${us}
fi

mv  plotLine_master_P_2_out_sub0.e plotLine_P2_${typem}_${ol}_$us.e
mv  plotLine_master_P_out_sub0.e plotLine_P1_${typem}_${ol}_$us.e
