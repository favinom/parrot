[Problem]
type = FEProblem
solve = false
[]

[Mesh]
 file = adapted_${resolution}_mesh.xda
 uniform_refine = ${unifSteps}
[]

[MeshModifiers]


[./aaassign1]
 type = AddSideSetsFromBoundingBox
 boundary_id_old = 'left'
 boundary_id_new = 13
 bottom_left = '-0.1  -0.1  -0.1'
 top_right =   '0.1  0.2500001 0.25000001'
 [../]
 
[./aaassign2]
 type = AddSideSetsFromBoundingBox
 boundary_id_old = 'bottom'
 boundary_id_new = 13
 bottom_left = '-0.1  -0.1  -0.1'
 top_right =   '0.250000001  0.1  0.250000001'
 [../]
 
[./aaassign3]
 type = AddSideSetsFromBoundingBox
 boundary_id_old = 'back'
 boundary_id_new = 13
 bottom_left = '-0.1  -0.1  -0.1'
 top_right =   '0.250000001  0.250000001  0.1'
 [../]

[./aaassign4]
type = AddSideSetsFromBoundingBox
boundary_id_old = 'right'
boundary_id_new = 10
bottom_left = '0 0.874999999 0.874999999'
top_right =   '1.00001  1.00001 1.00001'
[../]

[./aaassign5]
type = AddSideSetsFromBoundingBox
boundary_id_old = 'top'
boundary_id_new = 10
bottom_left = '0.874999999 0  0.874999999'
top_right =   '1.00001  1.00001 1.00001'
[../]

[./aaassign6]
type = AddSideSetsFromBoundingBox
boundary_id_old = 'front'
boundary_id_new = 10
bottom_left = '0.874999999  0.874999999 0'
top_right =   '1.00001  1.00001 1.00001'
[../]



[./aa1]
type = FractureUserObject
fn = 9
fx_string = '0.5,0.5,0.5,0.749975,0.75,0.749975,0.625,0.625,0.625'
fy_string = '0.5,0.5,0.5,0.749975,0.749975,0.75,0.625,0.625,0.625'
fz_string = '0.5,0.5,0.5,0.75,0.749975,0.749975, 0.625,0.625,0.625'
fa1_string = '0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0'
fa2_string = '0.0,90.0,0.0,0.0,90.0,0.0,0.0,90.0,0.0'
fa3_string = '0.0,0.0,90.0,0.0,0.0,90.0,0.0,0.0,90.0'
fd1_string = '1.0,1.0,1.0,0.50005,0.50005,0.50005,0.2501,0.2501,0.2501'
fd2_string = '1.0,1.0,1.0,0.50005,0.50005,0.50005,0.2501,0.2501,0.2501'
fd3_string = '0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001'
[../]

 [./assign1]
 type = SubdomainBoundingBox
 bottom_left = '0.50000001 -0.1 -0.1'
 top_right = '1.1 0.4999999 1.1'
 block_id = 1
 [../]
 
[./assign2]
 type = SubdomainBoundingBox
 block_id = 1
 bottom_left = '0.75000001 0.50000001 0.500000001'
 top_right = '1.1 0.749999999 1.1'
 [../]
 
 
[./assign3]
 type = SubdomainBoundingBox
 bottom_left = '0.62500001 0.500000001 0.5000000001'
 top_right = '0.749999999 0.624999999 0.7499999999'
 block_id = 1
 [../]


[./assign4]
type = RegionSubdomain
block_id = 2
regionMeshModifier = aa1
[../]
[]

[Variables]
[./pressure]  order = FIRST  family=LAGRANGE [../]
[]

#[AuxVariables]
#[./blockID] order=CONSTANT  family=MONOMIAL [../]
#[]

#[AuxKernels]
#[./ciao]
# type = MaterialRealAux
# variable = blockID
# property = RegionIDReal
# [../]
# []

#[Materials]
#[./hhhh]
#type = RegionMaterial
# regionMeshModifier = aa2
#[../]
#[]
 
[Kernels]
[./StressDivergenceParrot_real_x]
type = Reaction
variable = pressure
[../]
[]

[Executioner]
type=Steady
solve_type=LINEAR
line_search = 'none'
nl_abs_tol = 1e-8
[]

[Outputs]
file_base = refined_${resolution}_${unifSteps}
# exodus = true
print_linear_residuals = true
xdr = true
[]

[Adaptivity]
 marker = simplemark
 steps = ${adaptSteps}
 [./Markers]
 [./simplemark]
 type = BlockMarker
 blockNum = 2
 [../]
 [../]
 []
