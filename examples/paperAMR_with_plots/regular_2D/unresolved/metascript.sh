reslist="40 80" # "80 160 320 640"
aslist="8 9 10" # 4 8 16"

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
		./runScript.sh
	done
done
