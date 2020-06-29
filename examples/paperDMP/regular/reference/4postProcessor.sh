#!/bin/bash

source ./0defineVariable.sh

for f in DiffusionOutR_*.e;
do
	$pythonString ./scripts/postprocessor.py $f;
done

