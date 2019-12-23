typem=1;

ol=2;

np=4;

correction=1;

if [ $correction -eq 0 ]
then
	mpirun -n ${np} ../../../../parrot-opt -i 2advection.i typeMesh=${typem} origLevel=${ol}
fi
if [ $correction -eq 1 ]
then
	mpirun -n ${np} ../../../../parrot-opt -i 2advectionCorrection.i typeMesh=${typem} origLevel=${ol}
    mpirun -n ${np} ../../../../parrot-opt -i plotLine_master_P.i typeMesh=${typem} origLevel=${ol}
    mpirun -n ${np} ../../../../parrot-opt -i plotLine_master_P_2.i typeMesh=${typem} origLevel=${ol}
fi
