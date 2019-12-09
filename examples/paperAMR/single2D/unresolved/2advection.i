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
fn = 1
fx_string = '0.0'
fy_string = '0.0'
fa1_string = '-30.963756532073525'
fd1_string = '2000.0'
fd2_string = '0.01'
[../]
[]

[Variables]
[./CM] [../]
[]
 
[AuxVariables]
[./P_aux] [../]
[]

[Materials]
[./conductivity1] type = FractureMaterial fractureMeshModifier =  fractureUserObject block = 0
matrixPorosity = 0.2 fracturePorosity = 0.4
matrixPermeability = 1e-6 fracturePermeability = 1e-1
pressure = P_aux
[../]
[./conductivity2] type = FractureMaterial fractureMeshModifier =  fractureUserObject block = 2
matrixPorosity = 0.25 fracturePorosity = 0.4
matrixPermeability = 1e-5 fracturePermeability = 1e-1
pressure = P_aux
[../]
[]

[Kernels]
active='time upwind'
[upwind] type = Advection variable = CM [../]
[./time] type = PorosityTimeDerivative variable = CM lumping = true dim = 2 [../]
[]

[BCs]
[./u_injection_left] type = DirichletBC boundary = 11 variable = CM value='0.01' [../]
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
block_id='0 2'
value_p ='1e-6 1e-5 1e-1'
boundary_D_bc = '11 22'
value_D_bc='4.0 1.0'
boundary_N_bc = ''
value_N_bc=''
aux_variable=P_aux
fractureMeshModifier = fractureUserObject
#output_file=matrix.e
[../]
[]
