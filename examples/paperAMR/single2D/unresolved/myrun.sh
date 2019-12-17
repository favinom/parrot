#!/bin/bash
#!/bin/bash -l
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --constraint=mc
#SBATCH --job-name=test2
#SBATCH --time=0-24:00:00
#SBATCH --error=whale-opt_error.log
#SBATCH --account=u0
#SBATCH --partition=normal
#SBATCH --dependency=afterany:13036137
as=6;
np=4;
createmesh=1;
correction=1;

if [ $createmesh -eq 1 ]
then
	./myclean
	srun -n $SLURM_NTASKS --ntasks-per-node=12 -c $SLURM_CPUS_PER_TASK ../../../../parrot-opt -i 0refineBlock.i adapSteps=${as}
	for (( c=0; c<=as+1; c++ ))
	do
		if [ $c -le 9 ]
		then
			v=0$c
		else
			v=$c
		fi
		rm refinedMesh_00$v.xdr
	done

	for (( c=1; c<=as+1; c++ ))
	do
		temp=`expr $c - 1`
		if [ $c -le 9 ]
		then
			v1=0$c
		else
			v1=$c
		fi
		if [ $temp -le 9 ]
		then
			v2=0$temp
		else
			v2=$temp
		fi
		mv refinedMesh_00${v1}_mesh.xdr refinedMesh_00${v2}_mesh.xdr
	done
fi

if [ $as -le 9 ]
then
	as=0$as
else
	as=$as
fi
echo $as

#mpirun -n ${np} ../../../../parrot-opt -i 1diffusion.i adapSteps=${as}

if [ $correction -eq 0 ]
then
       srun -n $SLURM_NTASKS --ntasks-per-node=12 -c $SLURM_CPUS_PER_TASK ../../../../parrot-opt -i 2advection.i adapSteps=${as}
fi
if [ $correction -eq 1 ]
then
       srun -n $SLURM_NTASKS --ntasks-per-node=12 -c $SLURM_CPUS_PER_TASK ../../../../parrot-opt -i 2advectionCorrection.i adapSteps=${as}
fi
