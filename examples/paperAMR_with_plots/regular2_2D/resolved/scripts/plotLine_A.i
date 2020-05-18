[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1000
  xmin = -0.5
  xmax =  0.5
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
 vector_value = '0.5 0.5 0.0'
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
