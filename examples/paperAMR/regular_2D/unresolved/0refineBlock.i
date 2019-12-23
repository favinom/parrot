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
 nx = ${resolution}
 ny = ${resolution}
 dim = 2
[]

[MeshModifiers]
[./fractureUserObject]
 type = FractureUserObject
 fn = 6
 fx_string = '0.5,0.5,0.75,0.75,0.625,0.625'
 fy_string = '0.5,0.5,0.75,0.75,0.625,0.625'
 fa1_string = '0.0,90.0,0.0,90.0,0.0,90.0'
 fd1_string = '1.0,1.0,0.5,0.5,0.25,0.25'
 fd2_string = '1.0e-4,1.0e-4,1.0e-4,1.0e-4,1.0e-4,1.0e-4'
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
