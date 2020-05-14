declare -a adaptRefinement;
adaptRefinement[0]=" ' 14 0 ' "
#adaptRefinement[1]="'15 0'"
#adaptRefinement[2]="'16 0'"
#adaptRefinement[3]="'17 0'"
#adaptRefinement[4]="'18 0'"
declare -a alterRefinement;
alterRefinement[0]="'1 1 1 1 1 1 1 1 1 1 1 1 1 1'"         # 14
#alterRefinement[1]="'1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0'"     # 15
#alterRefinement[2]="'1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1'"     # 16
#alterRefinement[3]="'1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0'" # 17
#alterRefinement[4]="'1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1'" # 18
declare -a refinement;

refinementType=0
problemType=0
np=4

if [ $refinementType -eq 0 ]
then
	len=${#adaptRefinement[@]}
	for (( i=0; i<${len}; i++ ))
	do
		refinement[i]=${adaptRefinement[$i]}
		# echo ${refinement[$i]}
	done
fi

if [ $refinementType -eq 1 ]
then
	len=${#alterRefinement[@]}
	for (( i=0; i<${len}; i++ ))
	do
		refinement[i]=${alterRefinement[$i]}
		# echo ${refinement[$i]}
	done
fi

export np
export refinementType
export problemType

#./runScriptM.sh

len=${#refinement[@]}
for (( refLev=0; refLev<${len}; refLev++ )) # 
do
	refString=${refinement[$refLev]}
	export refString
	export refLev
	./postProcessorM.sh
done
