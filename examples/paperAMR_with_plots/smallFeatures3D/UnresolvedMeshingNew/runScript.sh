#!/bin/bash

correction=0;

if [ $as -le 9 ]
then
	as=0$as
else
	as=$as
fi
echo $as

if [ $correction -eq 0 ]
then
mpirun -n ${np} ../../../../parrot-opt -i  2advection.i  adapSteps=${as} 
fi

if [ $correction -eq 1 ]
then
mpirun -n ${np} ../../../../parrot-opt -i 2advectionCorrection.i adapSteps=${as}
fi


if [ $postprocessors -eq 1 ]
then
    mpirun -n ${np} ../../../../parrot-opt -i plotLine_master_P.i    adapSteps=${as}
    mpirun -n ${np} ../../../../parrot-opt -i plotLine_master_P_2.i  adapSteps=${as}
fi

mv plotLine_master_P_out_sub0.e    plotLine_p1_${as}.e
mv plotLine_master_P_2_out_sub0.e  plotLine_p2_${as}.e
