[Problem]
type = FEProblem
solve = false
[]

[Mesh]
type = GeneratedMesh
xmin= 0.0
xmax= 1.0
ymin= 0.0
ymax= 1.0
nx = ${mRes}
ny = ${mRes}
dim = 2
[]
 
 
[MeshModifiers]
 
[./subdomains_1]
 type = SubdomainBoundingBox
# top_right =   ' 0.5 1.001 0.0'
# bottom_left = '-0.01 -0.01 0.0'
 
 bottom_left = '-0.01 -0.01 0.0'
 top_right =   ' 1.01 0.7 0.0'
 
 
 block_id = 1
 block_name = block_1
 [../]
 
[./subdomains_2]
 type = SubdomainBoundingBox
# top_right =   '1.01 1.01 0.0'
# bottom_left = '0.5 -0.01 0.0'
 
 bottom_left = '0.0 0.7 0.0'
 top_right =   '1.01 1.01 0.0'
 
 
 block_id = 2
 block_name = block_2
 [../]
 
 
 
[./center_side_set_2]
 type = SideSetsBetweenSubdomains
 master_block = block_1
 paired_block = block_2
 new_boundary = 'new_side_set_2'
 depends_on = 'subdomains_1 subdomains_2'
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
file_base = scripts/refinedMesh
#exodus = true
perf_graph = true
xdr = true
[]

#[Adaptivity]
#marker = simplemark
#steps = ${adapSteps}
#[./Markers]
#[./simplemark]
#type = RegionMarker
#regionMeshModifier = fractureUserObject
#[../]
#[../]
#[]
