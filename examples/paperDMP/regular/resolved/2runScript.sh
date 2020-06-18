#!/bin/bash

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
	$clusterString $parrotString 2diffusion.i resolution=${res} unifSteps=${us} adaptSteps=${as}
fi
