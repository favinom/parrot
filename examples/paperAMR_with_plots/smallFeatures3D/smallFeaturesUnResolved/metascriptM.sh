#!/bin/bash -l
#BATCH --job-name=contact
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nestom@usi.ch
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --partition=fat
#SBATCH --mem=100000

back_x=18 #72
back_y=45  #160

aslist="3 4 5" # 3 4 5"

createmesh=1;
np=20;
postprocessors=1;

export back_x
export back_y
export np
as=5
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
