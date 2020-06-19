#!/bin/bash

# number of background elements
be=$1
# number of fracture elements
fe=$2

correction=$3

comman="runResolved"
source ./0defineVariable.sh

if [ $correction -eq 0 ]
then
	$clusterString $parrotString ./scripts/3advection.i mBe=${be} mFe=${fe}
fi
if [ $correction -eq 1 ]
then
	$clusterString $parrotString 3advection.i mBe=${be} mFe=${fe} UserObjects/active = 'soln storeOperatorsUO MassAssembly assembleVolumeVectors antidiffusiveFluxes'
fi
if [ $correction -eq 2 ]
then
	$clusterString $parrotString ./scripts/2diffusion.i mBe=${be} mFe=${fe}
fi
