problemType=0
np=1

declare -a amrList
declare -a umrList
resList[0]='10'; amrList[0]='1';   umrList[0]='0';
#resList[0]='100'; amrList[1]='10'   umrList[1]='0';
#resList[0]='100'; amrList[2]='11'   umrList[2]='0';
#resList[0]='100'; amrList[3]='10'   umrList[3]='1';

doRun=1;
doPost=1;

pythonString='/Applications/ParaView-5.5.1.app/Contents/bin/pvpython'
parrotString='mpirun -n '${np}' ../../../../parrot-opt -i '
clusterString=''
runString=$clusterString' '$parrotString
jobName=run_${amrName}_${umr}
errName=run_${amrName}_${umr}.err
outName=run_${amrName}_${umr}.out

function makeClusterStringLausanne
{
	runName=run_$1_$2_$3
	outName=out_$1_$2_$3.out
	errName=err_$1_$2_$3.err
	clusterString='bsub -q highmem -m node15 -n '${np}','${np}' -J '${runName}' -e '${errName}' -o '${outName}' '${np}' '${parrotString}
}  



len=${#amrList[@]}

for (( i=0; i<len; i++ )) # 
do
	export problemType
	export runString
	res=${resList[$i]}
	amr=${amrList[$i]}
	umr=${umrList[$i]}

	if [ $amr -le 9 ]
	then
		amrName=0$amr
	else
		amrName=$amr
	fi

	if [ $res -le 99 ]
	then
		resName=00$res
	else
		if [ $res -le 999 ]
		then
			resName=0$res
		else
			resName=$res
		fi
	fi

	export res
	export resName
	export amr
	export amrName
	export umr

	export runline

	if [ $doRun -eq 1 ]
	then
		./1runScript.sh
	fi
	if [ $doPost -eq 1 ]
	then
		 $pythonString scripts/postprocessor.py DiffusionOut2_${resName}_${amrName}_${umr}.e
	fi	
done


