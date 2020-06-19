#!/bin/bash

# number of background elements
be=$1
# number of fracture elements
fe=$2

#resolution
res=`expr ${be} + 3 \\* ${fe}`
echo $res

comman="generate"
source ./0defineVariable.sh
$clusterString $parrotString scripts/0generateMesh.i resolution=${res} mBe=${be} mFe=${fe}

comman="separate"
source ./0defineVariable.sh
$clusterString ./separateFile.sh uniform_${be}_${fe}.xda ${res}

comman="mat"
source ./0defineVariable.sh
commandString='faimesh('${be}','${fe}')'
$clusterString $matlabString ${commandString}

comman="unite"
source ./0defineVariable.sh
$clusterString ./uniteFile.sh adapted_${be}_${fe}.xda

comman="refine"
source ./0defineVariable.sh
$clusterString $parrotString scripts/1refineBlock.i mBe=${be} mFe=${fe}

