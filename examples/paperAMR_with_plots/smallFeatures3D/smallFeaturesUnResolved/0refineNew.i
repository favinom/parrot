[Problem]
type = FEProblem
solve = false
skip_nl_system_check = true
[]

[Mesh]
type = GeneratedMesh
dim = 3
xmin = 0.0
xmax = 1.0
ymin = 0.0
ymax = 2.25
zmin = 0.0
zmax = 1.0
nx = 9
ny = 20
nz = 9
parallel_type = DISTRIBUTED
[]

[MeshModifiers]
active = 'fractureUserObject'
[./createNewSidesetInflow]
type = AddSideSetsFromBoundingBox
boundary_id_old = 'bottom'
boundary_id_new = 11
top_right   = '10.1  1.9  0.666667'
bottom_left = '-0.1 -0.1  0.333332'
[../]

[./createNewSidesetOutflow1]
type = AddSideSetsFromBoundingBox
boundary_id_old = 'top'
boundary_id_new = 22
top_right =   ' 1.1 2.5 1.1'
bottom_left = '-0.1 1.8 0.665'
[../]

[./createNewSidesetOutflow2]
type = AddSideSetsFromBoundingBox
boundary_id_old = 'top'
boundary_id_new = 23
top_right =   ' 1.1 2.5 0.334'
bottom_left = '-0.1 1.8 0.0'
[../]

[./fractureUserObject]
type = FractureUserObject
fn = 8
fx_string = '0.5,0.5,0.5,0.5,0.2,0.2,0.77,0.83'
fy_string = '1.125,0.175,1.6,1.6,2.05,2.05,2.05,2.05'
fz_string = '0.5,0.5,0.675,0.31,0.5,0.5,0.5,0.5'
fd1_string = '0.9,0.9,0.9,0.9,0.4,0.4,0.4,0.4'
fd2_string = '1.75,0.25,1.25,1.2472,0.30594,0.30594,0.3,0.3'
fd3_string = '0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01'
fa1_string = '0,0,0,0,78.6901,-78.6901,0,0'
fa2_string = '0,90,0,0,-90,-90,-90,-90'
fa3_string = '0,0,16.2602,-15.8192,90,-90,0,0'
[../]

[]


[Executioner]
 type=Steady
 solve_type=LINEAR
 line_search = 'none'
 nl_abs_tol = 1e-8
 []

[UserObjects]
[./soln]
type = RefineMesh
fractureMeshModifier = fractureUserObject
refinements = '1 1 1 1 2 0'
filename = 'mesh.e'
[../]
[]