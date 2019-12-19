[Problem]
type = FEProblem
solve = false
[]

[Mesh]
 file = square_one_frac_${typem}.e
 block_id = '2 4 5 6 7'
 boundary_id = '1 2'
 boundary_name = 'inflow outflow'
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
 file_base = refinedMesh_${typem}_${unifSteps}
# exodus = true
 perf_graph = true
 xdr = true
[]

[Adaptivity]
 marker = simplemark
 steps = ${adapSteps}
 [./Markers]
 [./simplemark]
 type = BlockMarker
 blockNum = 4
 [../]
 [../]
[]
