createMesh=0
doRun=0
correction=0
doPost=1
np=4

declare -a beList
declare -a feList

beList[0]='160';   feList[0]='2';
#amrList[1]='10'   umrList[1]='0';
#amrList[2]='11'   umrList[2]='0';
#amrList[3]='10'   umrList[3]='1';

parrotString='mpirun -n ${np} ../../../../parrot-opt -i'

# matlabString='Users/mariagiuseppinanestola/Desktop/MATLAB_R2014b.app/bin/matlab -nodesktop -nosplash -r'
matlabString='/Applications/MATLAB_R2018b.app/bin/matlab -nodesktop -nosplash -r'
# matlabString='/soft/matlab/r2019a/bin/matlab -nodesktop -nosplash -r'

export parrotString
export matlabString

len=${#beList[@]}

for (( i=0; i<len; i++ )) # 
do
	be=${beList[$i]}
	fe=${feList[$i]}

	export be
	export fe

	jobName=run_${be}_${fe}
	errName=run_${be}_${fe}.err
	outName=run_${be}_${fe}.out

	# parrotString='bsub -q highmem -m node15 -n '${np}','${np}' -J '${jobName}' -e '${errName}' -o '${outName}' mpirun -n '${np}' ../../../../parrot-opt -i '
	# export parrotString

	if [ $createMesh -eq 1 ]
	then
		./1createMesh.sh
	fi

	if [ $doRun -eq 1 ]
	then
		export correction
		./2runScript.sh
	fi
	if [ $doPost -eq 1 ]
	then
		./3postProcessor.sh
	fi
#	if [ $doConv -eq 1 ]
#	then
#		./3conversionM.sh
#	fi
	
done

