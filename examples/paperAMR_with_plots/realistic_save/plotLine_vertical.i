[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1000
  xmin=-300
  xmax=300
[]
 
[MeshModifiers]
[./rotate]
 type = Transform
 transform = ROTATE
 vector_value = '0 0 90'
[../]

 [./translate]
 type = Transform
 transform = TRANSLATE
 vector_value = '625 300 0.0'
 [../]
[]


[Problem]
 type = FEProblem
 solve = false
 kernel_coverage_check = false
[]

 
[Variables]
  [./uv]
  [../]
[]
  
[Executioner]
 type = Transient
 num_steps = 1
 dt = 1
 
 solve_type = 'PJFNK'
 
# petsc_options_iname = '-pc_type -pc_hypre_type'
# petsc_options_value = 'hypre boomeramg'
 []

 
[Outputs]
[./out]
type = Exodus
execute_on = 'timestep_end'
[../]
[]
