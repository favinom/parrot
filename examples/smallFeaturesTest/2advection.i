[Problem]
type = ParrotProblem
use_AFC = true
[]

[Mesh]
 file = mesh_2_0.e
  boundary_id = '21 22'
  boundary_name = 'inflow outflow'
# partitioner = linear
[]

[MeshModifiers]
[./fractureUserObject]
type = FractureUserObject
fn = 8
fx_string = '0.5  ,0.5  ,0.77,0.83, 0.2,0.2  , 0.5,0.5'
fy_string = '1.125,0.175,2.05,2.05, 2.05,2.05 , 1.6,1.6'
fz_string = '0.5  ,0.5  ,0.5 ,0.5, 0.5,0.5 , 0.675,0.31'
fa1_string = '0,90,90,90,78.6901,-78.6901,0,0'
fa2_string = '0, 0, 0, 0,0,0,0,0'
fa3_string = '0,90,90,90,90,90,16.2602,-15.8192'
fd1_string = '0.9,0.25,0.3,0.3,0.3059,0.3059,0.9,0.9'
fd2_string = '1.75,0.9,0.4,0.4,0.4,0.4,1.25,1.2472'
fd3_string = '0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01'
[../]
[]

[Variables]
[./CM] [../]
[]
 
[AuxVariables]
[./P_aux] [../]
[]

[Materials]
[./conductivity1] type = FractureMaterial fractureMeshModifier =  fractureUserObject
matrixPorosity = 0.2 fracturePorosity = 0.2
matrixPermeability = 1 fracturePermeability = 1e4
pressure = P_aux
[../]
[]

[Kernels]
active='time upwind'
[upwind] type = Advection variable = CM [../]
[./time] type = PorosityTimeDerivative variable = CM lumping = true [../]
[]

[BCs]
[./u_injection_left] type = DirichletBC boundary = 21 variable = CM value='1' [../]
[]

[Preconditioning]
[./SMP]
type = SMP
full = true
[../]
[]

[Executioner]

type = Transient
solve_type= LINEAR
line_search = none

 petsc_options_iname=' -ksp_type            '   # -mat_view
 petsc_options_value='  ksp_parrot_preonly  '   # ::ascii_matlab

dt = 0.01
num_steps=100

[./Quadrature] order=TENTH [../]

[]

[Outputs]
 file_base = AdvectionOut_${adapSteps}
 exodus = true
 csv=true
 perf_graph = true
[]


[UserObjects]
[./soln]
type = SolveDiffusion
execute_on = 'initial'
block_id='1 2 3 4 5 6 7 11 12 13     '
value_p ='1 1 1 1 1 1 1 1  1  1  1e4 '
boundary_D_bc = '22'
value_D_bc='0.0'
boundary_N_bc = '21'
value_N_bc='-1.0'
aux_variable=P_aux
fractureMeshModifier = fractureUserObject
output_file=matrix.e
[../]
[]
