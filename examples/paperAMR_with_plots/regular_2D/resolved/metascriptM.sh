belist="80" # "80 160 320 640"
felist="1"

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
		./runScriptM.sh
	done
done
