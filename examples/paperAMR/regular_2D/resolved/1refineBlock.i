[Problem]
type = FEProblem
solve = false
[]

[Mesh]
 file = adapted_${resolution}_mesh.xda
 uniform_refine = ${unifSteps}
[]

 [MeshModifiers]
 
 [./aa1]
 type = FractureUserObject
 fn = 6
 fx_string = '0.5,0.5,0.75,0.75,0.625,0.625'
 fy_string = '0.5,0.5,0.75,0.75,0.625,0.625'
 fa1_string = '0.0,90.0,0.0,90.0,0.0,90.0'
 fd1_string = '1.0,1.0,0.5,0.5,0.25,0.25'
 fd2_string = '1.0e-4,1.0e-4,1.0e-4,1.0e-4,1.0e-4,1.0e-4'
 [../]

 [./aa2]
 type = BenchRegular2D
 [../]
 
 [./assign_id]
 type = RegionSubdomain
 block_id = 1
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
#exodus = true
print_linear_residuals = true
xdr = true
#xda=true
[]

[Adaptivity]
 marker = simplemark
 steps = ${adaptSteps}
 [./Markers]
 [./simplemark]
 type = BlockMarker
 blockNum = 1
 [../]
 [../]
 []
