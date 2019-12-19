[Problem]
type = FEProblem
solve = false
[]

[MeshModifiers]
 [./fractureUserObject]
 type = FractureUserObject
 fn = 2
 fx_string = '0.0,0.0'
 fy_string = '-6500.0,500.0'
 fa1_string = '79.695153531233970,-45.0'
 fd1_string = '20000000.0,200000000.0'
 fd2_string = '14.758048690536125,7.071067811865476'
 [../]
[]

 
[Mesh]
 file = hydrocoin.e
 uniform_refine = ${unifSteps}
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
 file_base = refinedMesh_${unifSteps}
# exodus = true
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
