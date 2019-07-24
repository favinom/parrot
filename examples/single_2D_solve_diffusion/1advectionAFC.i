[Problem]
type = ParrotProblem
use_AFC = true
[]

[Mesh]
 file = refinedMesh_${unifSteps}_000${adapSteps}_mesh.xdr
 block_id = '1 6 4 7 2'
 boundary_id = '6 7'
 boundary_name = 'inflow outflow'
 []

[Variables]
[./CM] [../]
[]
 
[AuxVariables]
[./pressure] [../]
[]

[Materials]
[./conductivity1] type = FlowAndTransport k = 1e-6 phi=0.2 block = '1 6 7' pressure = pressure [../]
[./conductivity4] type = FlowAndTransport k = 1e-1 phi=0.4 block = 4       pressure = pressure [../]
[./conductivity2] type = FlowAndTransport k = 1e-5 phi=0.25 block = 2      pressure = pressure [../]
 []

[Kernels]
active='time upwind'

[upwind]
type = Advection
variable = CM
[../]

[./time]
type = PorosityTimeDerivative
variable = CM
lumping = true
 dim = 2
[../]

[]

[BCs]
[./u_injection_left] type = DirichletBC boundary = inflow variable = CM value='0.01' [../]
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

 petsc_options_iname=' -ksp_type             '   # -mat_view
 petsc_options_value='  ksp_parrot_preonly    '   # ::ascii_matlab

dt = 1e7
num_steps=100

[./Quadrature]
order=TENTH
[../]

[]

[Outputs]
 file_base = AdvectionOut_${unifSteps}_${adapSteps}
exodus = true
csv=true
perf_graph = true
[]

[UserObjects]
[./myuo]
 type = SolveDiffusion
 execute_on = 'initial'
 block_id='1 2 4 6 7'
 value_p='1e-6 1e-5 1e-1 1e-6 1e-6'
 boundary_D_bc='6 7'
 function_D_bc='4.0 1.0'
 boundary_N_bc=' '
 value_N_bc=0.0
 [../]
 []
