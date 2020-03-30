
belist="80"
felist="1" # 4 8 16"


createmesh=0;

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
