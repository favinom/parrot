#!/bin/bash

#resolution
res=`expr ${be} + 3 \\* ${fe}`
echo $res
# number of adaptive steps

$parrotString scripts/0generateMesh.i resolution=${res} mBe=${be} mFe=${fe} --mesh-only
./separateFile.sh uniform_${be}_${fe}.xda ${res}

commandString='faimesh('${be}','${fe}')'

$matlabString ${commandString}

./uniteFile.sh adapted_${be}_${fe}.xda

$parrotString scripts/1refineBlock.i mBe=${be} mFe=${fe} --mesh-only

exit
# up to now, I just exit as we never used adaptivity for these simulations.
# We still leave this possibility

as=0;
us=0;

for (( c=0; c<=as+1; c++ ))
do
	rm refined_${res}_${us}_000$c.xdr
done

for (( c=1; c<=as+1; c++ ))
do
	temp=`expr $c - 1`
	mv refined_${res}_${us}_000${c}_mesh.xdr refined_${res}_${us}_000${temp}_mesh.xdr
done

for (( c=0; c<as; c++ ))
do
	rm refined_${res}_${us}_000${c}_mesh.xdr
done
