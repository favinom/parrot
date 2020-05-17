problemType=1
np=4

declare -a amrList
declare -a umrList
amrList[0]='9';   umrList[0]='0';
amrList[1]='10'   umrList[1]='0';
amrList[2]='11'   umrList[2]='0';
amrList[3]='10'   umrList[3]='1';

doRun=1;
doPost=0;
doConv=0;

len=${#amrList[@]}

for (( i=0; i<len; i++ )) # 
do
	export problemType
	amr=${amrList[$i]}
	umr=${umrList[$i]}

	if [ $amr -le 9 ]
	then
		amrName=0$amr
	else
		amrName=$amr
	fi

	export amr
	export amrName
	export umr

	exec='~/projects/parrot/parrot/parrot-opt'
	jobName=run_${amrName}_${umr}
	errName=run_${amrName}_${umr}.err
	outName=run_${amrName}_${umr}.out

	# # runline='bsub -q highmem -m node15 -n '${np}','${np}' -J '${jobName}' -e '${errName}' -o '${outName}' mpirun -n '${np}' ../../../parrot-opt -i '
	runline='mpirun -n '${np}' ../../../parrot-opt -i '
	export runline

	if [ $doRun -eq 1 ]
	then
		./1runScriptM.sh
	fi
	if [ $doPost -eq 1 ]
	then
		./2postProcessorM.sh
	fi
	if [ $doConv -eq 1 ]
	then
		./3conversionM.sh
	fi
	
done


