reslist="1" # 80 80 160" # 320 640"
aslist="15"

createmesh=1;

if [ $createmesh -eq 1 ]
then
	for res in $reslist
	do
		echo ${res}
		export res
		as=$aslist
		export as
		export createmesh
		./createMesh.sh
	done
fi

for res in $reslist
do
	for as in $aslist
	do
		export as
		export res
		./runScriptM.sh
	done
done
