problemType=0
np=4

doRun=1;
doPost=1;
doConv=1;

for (( refLev=5; refLev<=5; refLev++ )) # 
do
	
	export problemType

	if [ $refLev -le 9 ]
	then
		refLevName=0$refLev
	else
		refLevName=$refLev
	fi

	export refLev
	export refLevName

	exec='~/projects/parrot/parrot/parrot-opt'
	jobName=run_$refLevName
	errName=run_$refLevName.err
	outName=run_$refLevName.out

	# runline='bsub -q highmem -m node15 -n '${np}','${np}' -J '${jobName}' -e '${errName}' -o '${outName}' mpirun -n '${np}' ../../../parrot-opt -i '
	runline='mpirun -n '${np}' ../../../parrot-opt -i '
	export runline

	if [ $doRun -eq 1 ]
	then
		./1runScriptM.sh
	fi
	if [ $doRun -eq 1 ]
	then
		./2postProcessorM.sh
	fi
	if [ $doRun -eq 1 ]
	then
		./3conversionM.sh
	fi
	
done


