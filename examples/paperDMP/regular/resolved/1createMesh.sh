#!/bin/bash

# number of background elements
be=$1
# number of fracture elements
fe=$2

#resolution
res=`expr ${be} + 3 \\* ${fe}`
echo $res

source ./0defineVariable.sh

$clusterString $parrotString scripts/0generateMesh.i resolution=${res} mBe=${be} mFe=${fe}
$clusterString ./separateFile.sh uniform_${be}_${fe}.xda ${res}
commandString='faimesh('${be}','${fe}')'
$clusterString $matlabString ${commandString}
$clusterString ./uniteFile.sh adapted_${be}_${fe}.xda

$clusterString $parrotString scripts/1refineBlock.i mBe=${be} mFe=${fe}
