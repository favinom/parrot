[Problem]
type = FEProblem
solve = false
skip_nl_system_check = true
[]

[Mesh]
type = GeneratedMesh
dim = 2
xmin = 0.0
xmax = 2.0
ymin = 0.0
ymax = 2
nx = 1
ny = 1
parallel_type = distributed
[]


[MeshModifiers]

[./fractureUserObject]
type = FractureUserObject
fn = 2
fx_string = '1,1'
fy_string = '1,1'
fd1_string = '1,2'
fd2_string = '1,0.5'
fa1_string = '0,0'
[../]

[./my]
type = FractureRefinement
fractureMeshModifier = fractureUserObject
refinements = '1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0'
outputFileName= 'ciao.e'
doBoundaryRefinement = true
[../]

[]


[Executioner]
 type=Steady
 solve_type=LINEAR
 line_search = 'none'
 nl_abs_tol = 1e-8
[]
