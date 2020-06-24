[Mesh]
file = ../mesh_10_1.xdr
parallel_type = DISTRIBUTED
[]

[Problem]
solve = false
[]


[AuxVariables]
[./pressure_09_0] [../]
[./pressure_10_0] [../]
[./pressure_11_0] [../]
[./pressure_10_1] [../]
[]

[Functions]
active = ''
  [./pressureF_09_0]
    type = SolutionFunction
    scale_factor = 1.0
    solution = sol_09_0
  [../]

  [./pressureF_10_0]
    type = SolutionFunction
    scale_factor = 1.0
    solution = sol_10_0
  [../]

  [./pressureF_11_0]
    type = SolutionFunction
    scale_factor = 1.0
    solution = sol_11_0
  [../]

  [./pressureF_10_1]
    type = SolutionFunction
    scale_factor = 1.0
    solution = sol_10_1
  [../]

[]

 
[Preconditioning]
[./prec] type = SMP full = true ksp_norm = default [../]
[]
 
[Executioner]

 type=Steady
 solve_type= LINEAR
 line_search = none
# [./Quadrature] order = NINTH [../] # type = GRID
[]


[UserObjects]
[./sol_09_0]
type = SolutionUserObject
mesh = '../DiffusionOut_09_0.e'
system_variables = pressure
execute_on = 'initial'
timestep = 'LATEST'
[../]

[./sol_10_0]
type = SolutionUserObject
mesh = '../DiffusionOut_10_0.e'
system_variables = pressure
execute_on = 'initial'
timestep = 'LATEST'
[../]

[./sol_11_0]
type = SolutionUserObject
mesh = '../DiffusionOut_11_0.e'
system_variables = pressure
execute_on = 'initial'
timestep = 'LATEST'
[../]

[./sol_10_1]
type = SolutionUserObject
mesh = '../DiffusionOut_10_1.e'
system_variables = pressure
execute_on = 'initial'
timestep = 'LATEST'
[../]


[]

[AuxKernels]
[./sol_09_0]
type = SolutionAux
solution = sol_09_0
variable = pressure_09_0
scale_factor = 1.0
execute_on = 'initial'
[../]
[./sol_10_0]
type = SolutionAux
solution = sol_10_0
variable = pressure_10_0
scale_factor = 1.0
execute_on = 'initial'
[../]
[./sol_11_0]
type = SolutionAux
solution = sol_11_0
variable = pressure_11_0
scale_factor = 1.0
execute_on = 'initial'
[../]
[./sol_10_1]
type = SolutionAux
solution = sol_10_1
variable = pressure_10_1
scale_factor = 1.0
execute_on = 'initial'
[../]

[]



[Outputs]
 exodus = true
 file_base = differences
[]
