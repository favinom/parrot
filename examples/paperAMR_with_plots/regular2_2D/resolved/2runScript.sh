#!/bin/bash

if [ $correction -eq 0 ]
then
	$parrotString ./scripts/3advection.i mBe=${be} mFe=${fe}
fi
if [ $correction -eq 1 ]
then
	$parrotString 3advection.i mBe=${be} mFe=${fe} UserObjects/active = 'soln storeOperatorsUO MassAssembly assembleVolumeVectors antidiffusiveFluxes'
fi
if [ $correction -eq 2 ]
then
	$parrotString 2diffusion.i resolution=${res} unifSteps=${us} adaptSteps=${as}
fi
