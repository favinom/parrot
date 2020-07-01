#!/bin/bash
np=32
npReserve=32
node="-q highmem -m node15 "
#node="-m node05 "
jobName="${comman}_$1_$2"
outName="$jobName.out"
errName="$jobName.err"

echo "jobName= "$jobName

parrotString="mpirun -n ${np} ../../../../parrot-opt -i "

# matlabString='Users/mariagiuseppinanestola/Desktop/MATLAB_R2014b.app/bin/matlab -nodesktop -nosplash -r'
# matlabString="/Applications/MATLAB_R2018b.app/bin/matlab -nodesktop -nosplash -r "
matlabString='/soft/matlab/r2019a/bin/matlab -nodesktop -nosplash -r'

# pythonString="/Applications/ParaView-5.5.1.app/Contents/bin/pvpython "
pythonString="/home/mfavino/Applications/ParaView-5.8.0-MPI-Linux-Python3.7-64bit/bin/pvpython "

clusterString="bsub $node -n $npReserve,$npReserve -J ${jobName} -o ${outName} -e ${errName} "
# clusterString=" "

export parrotString
export matlabString
export pythonString
export clusterString
