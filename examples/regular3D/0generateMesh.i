[Problem]
 type = FEProblem
 solve = false
[]

[Mesh]
 type = GeneratedMesh
 dim = 3
 nx = ${resolution}
 ny = ${resolution}
 nz = ${resolution}
 xmin = 0
 xmax = 1
 ymin = 0
 ymax = 1
 zmin = 0
 zmax = 1
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

[Outputs]
 file_base =  uniform_${resolution}
 print_linear_residuals = false
 perf_graph = false
 xda = true
# exodus=true
[]
