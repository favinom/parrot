problemType=0
np=4

#declare -a amrList
#declare -a umrList
#declare -a resList

resList[0]='1600'; amrList[0]='6';   umrList[0]='0';
resList[1]='1600'; amrList[1]='7';   umrList[1]='0';
resList[2]='1600'; amrList[2]='8';   umrList[2]='0';
#
#
#resList[3]='400'; amrList[5]='6';   umrList[5]='0';
#resList[4]='400'; amrList[6]='7';   umrList[6]='0';
#resList[5]='400'; amrList[7]='8';   umrList[7]='0';
#
#
#resList[6]='800'; amrList[10]='6';  umrList[10]='0';
#resList[7]='800'; amrList[11]='7';  umrList[11]='0';
#resList[8]='800'; amrList[12]='8';  umrList[12]='0';
#
#resList[9]='1600'; amrList[10]='6';  umrList[10]='0';
#resList[10]='1600'; amrList[11]='7';  umrList[11]='0';
#resList[11]='1600'; amrList[12]='8';  umrList[12]='0';

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



len=${#amrList[0]}

#for (( i=0; i<len; i++ )) #
for (( i=0; i<3; i++ )) #
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

