typem=1;

ol=0; #2

np=4;

as=0;

us=0;

correction=0;

createmesh=0;

postprocessor=0;


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
	mpirun -n ${np} ../../../../parrot-opt -i 2advectionCorrection.i typeMesh=${typem} origLevel=${ol} adapSteps=${as} Uref=${us}
fi
if [ $postprocessor -eq 1 ]
then
    mpirun -n ${np} ../../../../parrot-opt -i plotLine_master_P.i typeMesh=${typem} origLevel=${ol} adapSteps=${as} Uref=${us}
    mpirun -n ${np} ../../../../parrot-opt -i plotLine_master_P_2.i typeMesh=${typem} origLevel=${ol} adapSteps=${as} Uref=${us}
fi
