np=4

pythonExec='/Applications/ParaView-5.5.1.app/Contents/bin/pvpython'

pythonScript='readlinesource_2D.py'

belist="1" # 80 160
AMRlist="12" # 9 10"

for be in $belist
do
	for AMR in $AMRlist
	do
		if [ $AMR -le 9 ]
		then
			v=0$AMR
		else
			v=$AMR
		fi

		#resolution
		# res=`expr ${be} + 3 \\* ${fe}`
		# echo $res
		mpirun -n ${np} ../../../parrot-opt -i plotLine_master_unresolved.i bepar=${be} amrpar=${v}
		$pythonExec $pythonScript horizontal_line_${be}_${v}.e
		$pythonExec $pythonScript vertical_line_${be}_${v}.e
	done
done
