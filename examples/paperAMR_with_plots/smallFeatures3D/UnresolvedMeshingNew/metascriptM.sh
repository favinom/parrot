#!/bin/bash -l

back_x=144  #288 #72 ##72
back_y=360  #720 #180  #160

aslist="2" # 3 4 5"

createmesh=1;
np=20;
postprocessors=0;

export back_x
export back_y
export np
as=2
export as

if [ $createmesh -eq 1 ]
then

	export createmesh
	./createMesh.sh
fi

for as in $aslist
do
	export postprocessors
	./runScript.sh
done