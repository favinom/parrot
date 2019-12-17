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
as=4;

as=5;
us=2;
np=1;
srun -n $SLURM_NTASKS --ntasks-per-node=12 -c $SLURM_CPUS_PER_TASK  ../../../parrot-opt -i 0refineBlock.i adapSteps=${as} unifSteps=${us}

for (( c=0; c<=as+1; c++ ))
do
rm refinedMesh_${us}_000$c.xdr
done

for (( c=1; c<=as+1; c++ ))
do
temp=`expr $c - 1`
mv refinedMesh_${us}_000${c}_mesh.xdr refinedMesh_${us}_000${temp}_mesh.xdr
done

srun -n $SLURM_NTASKS --ntasks-per-node=12 -c $SLURM_CPUS_PER_TASK  ../../../parrot-opt -i 1diffusion.i adapSteps=${as} unifSteps=${us}
#mpirun -n ${np} ../../parrot-opt -i 2advection.i adapSteps=${as} unifSteps=${us}
#mpirun -n ${np} ../../parrot-opt -i 2advectionAFC.i adapSteps=${as} unifSteps=${us}
