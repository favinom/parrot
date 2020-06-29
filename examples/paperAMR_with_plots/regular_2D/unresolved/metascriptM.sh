reslist="80" # 80 80 160" # 320 640"
aslist="7"

createmesh=0;

if [ $createmesh -eq 1 ]
then
	for res in $reslist
	do
		echo ${res}
		export res
		as=10
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
