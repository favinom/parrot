#!/bin/bash
np=16
npReserve=16
nod="-q highmem -m node15 "
jobName="run_$1_$2"
outName="$jobName.out"
errName="$jobName.err"

parrotString="mpirun -n ${np} ../../../../parrot-opt -i "

# matlabString='Users/mariagiuseppinanestola/Desktop/MATLAB_R2014b.app/bin/matlab -nodesktop -nosplash -r'
matlabString="/Applications/MATLAB_R2018b.app/bin/matlab -nodesktop -nosplash -r "
# matlabString='/soft/matlab/r2019a/bin/matlab -nodesktop -nosplash -r'

pythonString="/Applications/ParaView-5.5.1.app/Contents/bin/pvpython "
# pythonString='~/Applications/ParaView-5.8.0-MPI-Linux-Python3.7-64bit/bin/pvpython'

# clusterString="bsub $node -n $npReserve,$npReserve -J ${jobName} -o ${outName} -e ${errName} "
clusterString=" "

export parrotString
export matlabString
export pythonString
export clusterString
