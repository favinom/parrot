problemType=0
np=16

declare -a amrList
declare -a umrList
declare -a resList
resList[0]='160'; amrList[0]='6';   umrList[0]='0';
resList[1]='160'; amrList[1]='7';   umrList[1]='0';
resList[2]='160'; amrList[2]='8';   umrList[2]='0';
resList[3]='160'; amrList[3]='9';   umrList[3]='0';

resList[4]='320'; amrList[4]='5';   umrList[4]='0';
resList[5]='320'; amrList[5]='6';   umrList[5]='0';
resList[6]='320'; amrList[6]='7';   umrList[6]='0';
resList[7]='320'; amrList[7]='8';   umrList[7]='0';

resList[8]='640'; amrList[8]='4';   umrList[8]='0';
resList[9]='640'; amrList[9]='5';   umrList[9]='0';
resList[10]='640'; amrList[10]='6';   umrList[10]='0';
resList[11]='640'; amrList[11]='7';   umrList[11]='0';

resList[12]='1280'; amrList[12]='3';   umrList[12]='0';
resList[13]='1280'; amrList[13]='4';   umrList[13]='0';
resList[14]='1280'; amrList[14]='5';   umrList[14]='0';
resList[15]='1280'; amrList[15]='6';   umrList[15]='0';

doRun=1;
doPost=0;

#pythonString='/Applications/ParaView-5.5.1.app/Contents/bin/pvpython'
pythonString='~/Applications/ParaView-5.8.0-MPI-Linux-Python3.7-64bit/bin/pvpython'
parrotString='mpirun -n '${np}' ../../../../parrot-opt -i '
clusterString=''

function makeClusterStringLausanne
{
	runName=run_$1_$2_$3
	outName=out_$1_$2_$3.out
	errName=err_$1_$2_$3.err
	# clusterString='bsub -q highmem -m node15 -n '$4','$4' -J '${runName}' -e '${errName}' -o '${outName}' '${np}' '
	clusterString='bsub -m node07 -n '$4','$4' -J '${runName}' -e '${errName}' -o '${outName}' '
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

	makeClusterStringLausanne $res $amr $umr 16	
        runString=$clusterString' '$parrotString
	export runString
	echo $runString
	
	if [ $doRun -eq 1 ]
	then
		./1runScript.sh
	fi
	if [ $doPost -eq 1 ]
	then
		 $pythonString scripts/postprocessor.py DiffusionOutN_${resName}_${amrName}_${umr}.e
	fi	
done


