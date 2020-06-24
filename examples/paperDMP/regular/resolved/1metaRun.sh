#!/bin/bash

# action="./2createMesh.sh"
action="./3runScript.sh"

correction=2

$action 5120 64 $correction
exit

if [ 0 -gt 100 ]
then

$action   80  1 $correction
$action   80  2 $correction
$action   80  4 $correction
$action   80  8 $correction
$action   80 16 $correction

$action  160  1 $correction
$action  160  2 $correction
$action  160  4 $correction
$action  160  8 $correction
$action  160 16 $correction

$action  320  1 $correction
$action  320  2 $correction
$action  320  4 $correction
$action  320  8 $correction
$action  320 16 $correction

$action  640  1 $correction
$action  640  2 $correction
$action  640  4 $correction
$action  640  8 $correction
$action  640 16 $correction

$action 1280  1 $correction
$action 1280  2 $correction
$action 1280  4 $correction
$action 1280  8 $correction
$action 1280 16 $correction

fi

$action 2560  1 $correction
$action 2560  2 $correction
$action 2560  4 $correction
$action 2560  8 $correction
$action 2560 16 $correction

$action 5120  1 $correction
$action 5120  2 $correction
$action 5120  4 $correction
$action 5120  8 $correction
$action 5120 16 $correction

exit
$action 4000  1 $correction
$action 4000  2 $correction
$action 4000  4 $correction
$action 4000  8 $correction
$action 4000 16 $correction

exit

$action 8000  1 $correction
$action 8000  2 $correction
$action 8000  4 $correction
$action 8000  8 $correction
$action 8000 16 $correction
