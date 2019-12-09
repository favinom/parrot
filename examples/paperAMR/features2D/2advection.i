[Problem]
type = ParrotProblem
use_AFC = true
[]

[Mesh]
 file = refinedMesh_00${adapSteps}_mesh.xdr
 boundary_id = '11 22'
 boundary_name = 'inflow outflow'
[]

[MeshModifiers]
[./fractureUserObject]
type = FractureUserObject
fn = 7
fx_string = '0.0,50,-60,-50,-100.0,-125.0,-150.0'
fy_string = '100.0,-100,-100,0.0,25.0,50.0,-75.0'
fa1_string = '90.0,-63.4349488229,-116.5650511771,0.0,90.0,0.0,90.0'
fd1_string = '200.0,223.60679775,223.60679775,100.0,50.0,50.0,250.0'
fd2_string = '0.5,0.5,0.5,0.5,0.5,0.5,0.5'
[../]
[]

[Variables]
[./CM] [../]
[]
 
[AuxVariables]
[./P_aux] [../]
[]

[Materials]
[./conductivity1] type = FractureMaterial fractureMeshModifier =  fractureUserObject matrixPorosity = 0.2 fracturePorosity = 0.4
matrixPermeability = 1e-6 fracturePermeability = 1e-1
pressure = P_aux
[../]
[]

[Kernels]
active='time upwind'
[upwind] type = Advection variable = CM [../]
[./time] type = PorosityTimeDerivative variable = CM lumping = true [../]
[]

[BCs]
[./u_injection_left] type = DirichletBC boundary = top variable = CM value='0.01' [../]
#[./u_injection_left] type = NeumannBC boundary = top variable = CM value='1.0' [../]
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

dt = 1e7
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
block_id='0'
value_p ='1e-6 1e-1'
boundary_D_bc = ' 0    2   '
value_D_bc=     ' 1.0  4.0 '
boundary_N_bc = ''
value_N_bc=''
aux_variable=P_aux
fractureMeshModifier = fractureUserObject
#output_file=matrix.e
[../]
[]
