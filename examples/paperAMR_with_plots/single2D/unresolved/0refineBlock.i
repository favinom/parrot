[Problem]
type = FEProblem
solve = false
[]

[Mesh]
type = GeneratedMesh
xmin=-50.0
xmax= 50.0
ymin=-50.0
ymax= 50.0
nx = 1600
ny = 1600
dim = 2
[]

[MeshModifiers]
[./createNewSidesetInflow]
type = AddSideSetsFromBoundingBox
boundary_id_old = 'left'
boundary_id_new = 11
top_right   = '-45.0  51.0  0.0'
bottom_left = '-51.0  40.0  0.0'
[../]

[./createNewSidesetOutflow]
type = AddSideSetsFromBoundingBox
boundary_id_old = 'right'
boundary_id_new = 22
bottom_left = '45.0 -51.0 0.0'
top_right =   '51.0 -40.0 0.0'
[../]
 
[./subdomains_0]
type = SubdomainBoundingBox
bottom_left = '-51.00 -40.00001 0.0'
top_right =   '0.0000 51.0001 0.0'
block_id = 0
block_name = up_block_0
[../]

[./subdomains_1]
type = SubdomainBoundingBox
bottom_left = '-51.00 -51.0 0.0'
top_right =   ' 0.00 -40.0 0.0'
block_id = 2
block_name = bottom_block_1
[../]

[./subdomains_2]
type = SubdomainBoundingBox
bottom_left = '0.00001 -51.0 0.0'
top_right =   '51.00000 -40.0 0.0'
block_id = 1
block_name = bottom_block_2
[../]


[./subdomains_3]
type = SubdomainBoundingBox
bottom_left = '0.00001 -40.00001 0.0'
top_right =   '51.000 51.0001 0.0'
block_id = 3
block_name = up_block_3
[../]

[./center_side_set_1]
type = SideSetsBetweenSubdomains
master_block = up_block_0
paired_block = up_block_3
new_boundary = 'new_side_set_1'
depends_on = 'subdomains_3 subdomains_0'
[../]

[./center_side_set_2]
type = SideSetsBetweenSubdomains
master_block = bottom_block_1
paired_block = bottom_block_2
new_boundary = 'new_side_set_2'
depends_on = 'subdomains_1 subdomains_2'
[../]
 
[./fractureUserObject]
type = FractureUserObject
fn = 1
fx_string = '0.0'
fy_string = '0.0'
fa1_string = '-30.963756532073525'
fd1_string = '2000.0'
fd2_string = '0.01'
[../]

[]

[Variables]
[./pressure]  order = FIRST  family=LAGRANGE [../]
[]

[Kernels]
[./StressDivergenceParrot_real_x] type = Reaction variable = pressure [../]
[]

[Executioner]
type=Steady
solve_type=LINEAR
line_search = 'none'
nl_abs_tol = 1e-8
[]

[Outputs]
file_base = refinedMesh
#exodus = true
perf_graph = true
xdr = true
[]

[Adaptivity]
marker = simplemark
steps = ${adapSteps}
[./Markers]
[./simplemark]
type = RegionMarker
regionMeshModifier = fractureUserObject
[../]
[../]
[]
