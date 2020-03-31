[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1000
  xmin=-0.5
  xmax=0.5
[]
 
[MeshModifiers]
[./rotate]
 type = Transform
 transform = ROTATE
 vector_value = '90 90 0'
[../]

 [./translate]
 type = Transform
 transform = TRANSLATE
 vector_value = '0.5 1.1 0.5'
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
[./out]
 execute_on = 'timestep_end'
 type = Exodus
[../]
[]
