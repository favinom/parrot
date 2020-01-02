typem=1;

ol=2;

np=4;

as=0;

correction=0;

createmesh=1;
postprocessor=1;


if [ $createmesh -eq 1 ]
then
 mpirun -n ${np} ../../../../parrot-opt -i 0refineBlock.i typeMesh=${typem} origLevel=${ol} adapSteps=${as}
fi

if [ $correction -eq 0 ]
then
	mpirun -n ${np} ../../../../parrot-opt -i 2advection.i typeMesh=${typem} origLevel=${ol} adapSteps=${as}
fi
if [ $correction -eq 1 ]
then
	mpirun -n ${np} ../../../../parrot-opt -i 2advectionCorrection.i typeMesh=${typem} origLevel=${ol} adapSteps=${as}
fi
if [ $postprocessor -eq 1 ]
then
    mpirun -n ${np} ../../../../parrot-opt -i plotLine_master_P.i typeMesh=${typem} origLevel=${ol} adapSteps=${as}
    mpirun -n ${np} ../../../../parrot-opt -i plotLine_master_P_2.i typeMesh=${typem} origLevel=${ol} adapSteps=${as}
fi
