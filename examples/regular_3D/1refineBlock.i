[Problem]
type = FEProblem
solve = false
[]

[Mesh]
 file = adapted_${resolution}_mesh.xda
 uniform_refine = ${unifSteps}
[]

 [MeshModifiers]
# active=''
 [./bb_omega31]  type = Omega3Regular3D [../]

 [./bb_fractures]
 type = FractureUserObject
 fn = 9
 fx_string = '0.5,0.5,0.5,0.749975,0.749975,0.749975,0.625,0.625,0.625'
 fy_string = '0.5,0.5,0.5,0.749975,0.749975,0.749975,0.625,0.625,0.625'
 fz_string = '0.5,0.5,0.5,0.749975,0.749975,0.749975,0.625,0.625,0.625'
 fa1_string = '0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0'
 fa2_string = '0.0,90.0,0.0,0.0,90.0,0.0,0.0,90.0,0.0'
 fa3_string = '0.0,0.0,90.0,0.0,0.0,90.0,0.0,0.0,90.0'
 fd1_string = '1.0,1.0,1.0,0.50005,0.50005,0.50005,0.2501,0.2501,0.2501'
 fd2_string = '1.0,1.0,1.0,0.50005,0.50005,0.50005,0.2501,0.2501,0.2501'
 fd3_string = '1.0e-4,1.0e-4,1.0e-4,1.0e-4,1.0e-4,1.0e-4,1.0e-4,1.0e-4,1.0e-4'
 [../]

 [./aa_assign_id]
 type = RegionSubdomain
 block_id = 2
 regionMeshModifier = bb_omega31
 [../]

 [./bb_assign_id]
 type = RegionSubdomain
 block_id = 1
 regionMeshModifier = bb_fractures
 [../]

 [./aaa_par_omega_f_0]
 type = AddSideSetsFromBoundingBox
 boundary_id_old = 'left'
 boundary_id_new = 13
 bottom_left = '-0.1  -0.1  -0.1'
 top_right =   '0.1  0.2500001 0.25000001'
 [../]
 
[./aaa_par_omega_f_1]
 type = AddSideSetsFromBoundingBox
 boundary_id_old = 'bottom'
 boundary_id_new = 13
 bottom_left = '-0.1  -0.1  -0.1'
 top_right =   '0.250000001  0.1  0.250000001'
 [../]
 
[./aaa_par_omega_f_2]
 type = AddSideSetsFromBoundingBox
 boundary_id_old = 'back'
 boundary_id_new = 13
 bottom_left = '-0.1  -0.1  -0.1'
 top_right =   '0.250000001  0.250000001  0.1'
 [../]
 
[./aaa_par_omega_h_0]
 type = AddSideSetsFromBoundingBox
 boundary_id_old = 'right'
 boundary_id_new = 14
 bottom_left = '0 0.874999999 0.874999999'
 top_right =   '1.00001  1.00001 1.00001'
 [../]
 
[./aaa_par_omega_h_1]
 type = AddSideSetsFromBoundingBox
 boundary_id_old = 'top'
 boundary_id_new = 14
 bottom_left = '0.874999999 0  0.874999999'
 top_right =   '1.00001  1.00001 1.00001'
 [../]
 
[./aaa_par_omega_h_2]
 type = AddSideSetsFromBoundingBox
 boundary_id_old = 'front'
 boundary_id_new = 14
 bottom_left = '0.874999999  0.874999999 0'
 top_right =   '1.00001  1.00001 1.00001'
 [../]

 []

[Variables]
[./pressure]  order = FIRST  family=LAGRANGE [../]
[]

[AuxVariables]
[./howMany] order=CONSTANT  family=MONOMIAL [../]
[]

[AuxKernels]
[./ciao]
 type = MaterialRealAux
 variable = howMany
 property = HowManyFractures
 [../]
 []

[Materials]
[./hhhh] type = RegionMaterial regionMeshModifier = bb_fractures countFractures = true [../]
[]
 
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
 [./Quadrature] order = constant [../] # 
[]

[Outputs]
 file_base = refined_${resolution}_${unifSteps}
exodus = true
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
