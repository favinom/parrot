
#!/bin/bash


if [ $createmesh -eq 1 ]
then
#	./myclean.sh
	mpirun -n ${np}  ../../../../parrot-opt -i 0refineBlock.i adapSteps=${as} nx=${back_x} ny=${back_y}
	for (( c=0; c<=as+1; c++ ))
	do
		if [ $c -le 9 ]
		then
			v=0$c
		else
			v=$c
		fi
		rm refinedMesh_00$v.xdr
	done


         mv refinedMesh_0000_mesh.xdr initial_mesh.xdr

fi
