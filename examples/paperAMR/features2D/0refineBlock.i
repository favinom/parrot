[Problem]
type = FEProblem
solve = false
[]

[Mesh]
 type = GeneratedMesh
 xmin=-200.0
 xmax= 200.0
 ymin=-200.0
 ymax= 200.0
 nx = 32
 ny = 32
 dim = 2
[]

[MeshModifiers]
[./fractureUserObject]
type = FractureUserObject
fn = 7
fx_string = '0.0,50,-60,-50,-100.0,-125.0,-150.0'
fy_string = '100.0,-100,-100,0.0,25.0,50.0,-75.0'
fa1_string = '90.0,-63.4349488229,-116.5650511771,0.0,90.0,0.0,90.0'
fd1_string = '200.0,223.60679775,223.60679775,100.0,50.0,50.0,250.0'
fd2_string = '0.5,0.5,0.5,0.5,0.5,0.5,0.5'
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
