[Problem]
type = FEProblem
solve = false
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
[]

[MeshModifiers]
#active =''
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
[]

[Executioner]
 type=Steady
 solve_type=LINEAR
 line_search = 'none'
 nl_abs_tol = 1e-8
 []


[Outputs]
 file_base = refinedMesh
 xdr=true
[]
