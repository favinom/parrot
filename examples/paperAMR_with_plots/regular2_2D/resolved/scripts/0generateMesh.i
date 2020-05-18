[Problem]
 type = FEProblem
 solve = false
[]

[Mesh]
 type = GeneratedMesh
 dim = 2
 nx = ${resolution}
 ny = ${resolution}
 xmin = 0
 xmax = 1
 ymin = 0
 ymax = 1
# partitioner = parmetis
 parallel_type = distributed
[]

[MeshModifiers]

[./fractureUserObject]
type = FractureUserObject
fn = 1
fx_string  = '0.0'
fy_string  = '0.0'
fa1_string = '0.0'
fd1_string = '0.0'
fd2_string = '0.0'
[../]


[./my]
type = FractureRefinement
fractureMeshModifier = fractureUserObject
refinements='0 0' 
outputFileName = uniform_${mBe}_${mFe}.xda
# doBoundaryRefinement = true
[../]
[]

[Executioner]
 type=Steady
 solve_type=LINEAR
 line_search = 'none'
[]

