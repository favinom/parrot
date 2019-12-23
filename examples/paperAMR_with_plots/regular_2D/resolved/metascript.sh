belist="80 160 320 640"
felist="1 2 4 8" # 4 8 16"

createmesh=1;

for be in $belist
do
	for fe in $felist
	do
		export fe
		export be
		if [ $createmesh -eq 1 ]
		then
			./createmesh.sh
		fi
		./runScript.sh
	done
done
