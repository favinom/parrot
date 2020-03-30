[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1000
  xmin=0
  xmax=1
[]
 
[MeshModifiers]
 [./translate]
 type = Transform
 transform = TRANSLATE
 vector_value = '0 2.15 0.5'
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
 type = Exodus
[]
