[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1000
  xmin=-70.7106781187
  xmax=70.7106781187
[]
 
[MeshModifiers]
 active='rotate1'
 [./rotate1]
 type = Transform
 transform = ROTATE
 vector_value = '-30.963756532073525 0 0'
[../]
 [./rotate2]
 type = Transform
 transform = ROTATE
 vector_value = '-45 0 0'
 [../]
 []


[Problem]
 type = FEProblem
 solve = false
 kernel_coverage_check = false
[]

 
[Variables]
  [./u]
  [../]
[]

[AuxVariables]
  [./v]
  [../]
[]
 
 
[Executioner]
 type = Transient
 num_steps = 1
 dt = 1
 
 solve_type = 'PJFNK'
 
 petsc_options_iname = '-pc_type -pc_hypre_type'
 petsc_options_value = 'hypre boomeramg'
 []

 
[Outputs]
 execute_on = 'timestep_end'
 [./out]
 type = Exodus
 execute_on = 'timestep_end'
 [../]
[]
