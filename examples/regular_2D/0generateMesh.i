
[Mesh]
 type = GeneratedMesh
 dim = 2
 nx = ${resolution}
 ny = ${resolution}
 xmin = 0
 xmax = 1
 ymin = 0
 ymax = 1
# elem_type = QUAD9
# partitioner = parmetis
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
[]

[Problem]
 type = FEProblem
 solve = false
[]

[Outputs]
 file_base =  uniform_${resolution}
 print_linear_residuals = false
 perf_graph = false
 xda = true
# exodus=true
[]
