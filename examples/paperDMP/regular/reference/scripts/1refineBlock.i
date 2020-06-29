[Problem]
type = FEProblem
solve = false
[]

[Variables]
[./dummy] [../]
[]

[Mesh]
 file = ../adapted_${mBe}_${mFe}.xda
 parallel_type = distributed
[]


[MeshModifiers]
[./subdomains_1]
type = SubdomainBoundingBox
bottom_left = '-0.01 -0.01 0.0'
top_right =   ' 1.01 0.7 0.0'
block_id = 1
block_name = block_1
[../]

[./subdomains_2]
type = SubdomainBoundingBox
bottom_left = '0.0 0.7 0.0'
top_right =   '1.01 1.01 0.0'
block_id = 2
block_name = block_2
[../]

[./zzzcenter_side_set_2]
type = SideSetsBetweenSubdomains
master_block = block_1
paired_block = block_2
new_boundary = 'new_side_set_2'
depends_on = 'subdomains_1 subdomains_2'
[../]
 
[./fractureUserObject]
type = FractureUserObject
fn = 6
fx_string = '0.5,0.5,0.75,0.75,0.625,0.625'
fy_string = '0.5,0.5,0.75,0.75,0.625,0.625'
fa1_string = '0.0,90.0,0.0,90.0,0.0,90.0'
fd1_string = '1.0,1.0,0.5,0.5,0.25,0.25'
fd2_string = '1.0e-4,1.0e-4,1.0e-4,1.0e-4,1.0e-4,1.0e-4'
[../]

[./output]
type = FractureRefinement
fractureMeshModifier = fractureUserObject
refinements='0 0' 
outputFileName = refined_${mBe}_${mFe}.xdr
# doBoundaryRefinement = true
depends_on = 'zzzcenter_side_set_2'
[../]

[]

[Executioner]
type=Steady
solve_type=LINEAR
[]
