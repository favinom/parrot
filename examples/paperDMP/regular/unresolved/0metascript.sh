problemType=1
np=4

declare -a amrList
declare -a umrList
declare -a resList
resList[0]='80'; amrList[0]='7';   umrList[0]='0';
resList[1]='80'; amrList[1]='8';   umrList[1]='0';
resList[2]='80'; amrList[2]='9';   umrList[2]='0';
resList[3]='80'; amrList[3]='10';   umrList[3]='0';
resList[4]='80'; amrList[4]='11';   umrList[4]='0';

resList[5]='160'; amrList[5]='6';   umrList[5]='0';
resList[6]='160'; amrList[6]='7';   umrList[6]='0';
resList[7]='160'; amrList[7]='8';   umrList[7]='0';
resList[8]='160'; amrList[8]='9';   umrList[8]='0';
resList[9]='160'; amrList[9]='10';  umrList[9]='0';

resList[10]='320'; amrList[10]='5';  umrList[10]='0';
resList[11]='320'; amrList[11]='6';  umrList[11]='0';
resList[12]='320'; amrList[12]='7';  umrList[12]='0';
resList[13]='320'; amrList[13]='8';  umrList[13]='0';
resList[14]='320'; amrList[14]='9';  umrList[14]='0';

resList[15]='640'; amrList[15]='4';  umrList[15]='0';
resList[16]='640'; amrList[16]='5';  umrList[16]='0';
resList[17]='640'; amrList[17]='6';  umrList[17]='0';
resList[18]='640'; amrList[18]='7';  umrList[18]='0';
resList[19]='640'; amrList[19]='8';  umrList[19]='0';

resList[20]='1280'; amrList[20]='3';   umrList[20]='0';
resList[21]='1280'; amrList[21]='4';   umrList[21]='0';
resList[22]='1280'; amrList[22]='5';   umrList[22]='0';
resList[23]='1280'; amrList[23]='6';   umrList[23]='0';
resList[24]='1280'; amrList[24]='7';   umrList[24]='0';

resList[25]='2560'; amrList[25]='2';   umrList[25]='0';
resList[26]='2560'; amrList[26]='3';   umrList[26]='0';
resList[27]='2560'; amrList[27]='4';   umrList[27]='0';
resList[28]='2560'; amrList[28]='5';   umrList[28]='0';
resList[29]='2560'; amrList[29]='6';   umrList[29]='0';

doRun=1;
doPost=0;

#pythonString='/Applications/ParaView-5.5.1.app/Contents/bin/pvpython'
pythonString='~/Applications/ParaView-5.8.0-MPI-Linux-Python3.7-64bit/bin/pvpython '
parrotString='mpirun -n '${np}' ../../../../parrot-opt -i '
clusterString=''

function makeClusterStringLausanne
{
	runName=runRS_$1_$2_$3
	outName=runRS_$1_$2_$3.out
	errName=runRS_$1_$2_$3.err
	clusterString='bsub -q highmem -m node15 -n '$4','$4' -J '${runName}' -e '${errName}' -o '${outName}' '
	# clusterString='bsub -m node10 -n '$4','$4' -J '${runName}' -e '${errName}' -o '${outName}' '
}  



len=${#amrList[@]}

#for (( i=0; i<len; i++ )) # 
for (( i=2; i<3; i++ )) # 
do
	export problemType
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

	#makeClusterStringLausanne $res $amr $umr 16
        runString=$clusterString' '$parrotString
	export runString
	if [ $doRun -eq 1 ]
	then
		echo $runString
		./1runScript.sh
	fi
	if [ $doPost -eq 1 ]
	then
		makeClusterStringLausanne $res $amr $umr 1
        	$clusterString $pythonString scripts/postprocessor.py DiffusionOutN_${resName}_${amrName}_${umr}.e
	fi	
done

