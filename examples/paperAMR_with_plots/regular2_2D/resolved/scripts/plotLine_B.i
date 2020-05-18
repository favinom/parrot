[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1000
  xmin = 0
  xmax = 1
[]
 
[MeshModifiers]

 [./translate]
 type = Transform
 transform = TRANSLATE
 vector_value = '0.0 0.75 0.0'
 [../]
[]


[Problem]
 type = FEProblem
 solve = false
 kernel_coverage_check = false
[]

 
[Variables]
[./p] [../]
[./cm1] [../]
[./cm2] [../]
[./cm3] [../]
[]
  
[Executioner]
 type = Transient
 num_steps = 1
 dt = 1
 
 []

 
[Outputs]
[./out]
type = Exodus
execute_on = 'timestep_end'
[../]
[]
