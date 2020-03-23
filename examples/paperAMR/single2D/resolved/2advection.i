[Problem]
type = ParrotProblem
use_AFC = true
operator_userobject = storeOperatorsUO
[]

[GlobalParams]
conservative = true
[]


[Mesh]
file = refinedMesh_${unifSteps}_000${adapSteps}_mesh.xdr
block_id = '1 6 4 7 2'
boundary_id = '6 7'
boundary_name = 'inflow outflow'
uniform_refine = 0
[]

[Variables]
[./CM] [../]
[]

[AuxVariables]
[./P_aux] [../]
[]

[Materials]
[./conductivity1] type = FlowAndTransport k = 1e-6 phi=0.2  conservative = true block = '1 6 7' pressure = P_aux [../]
[./conductivity4] type = FlowAndTransport k = 1e-1 phi=0.4  conservative = true block = 4       pressure = P_aux [../]
[./conductivity2] type = FlowAndTransport k = 1e-5 phi=0.25 conservative = true block = 2       pressure = P_aux [../]
[]


[Kernels]
active='time upwind'
[upwind] type = Advection variable = CM [../]
[./time] type = PorosityTimeDerivative variable = CM lumping = true dim = 2 [../]
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

[./TimeIntegrator]
type = VoidIntegrator
[../]


type = Transient
solve_type= LINEAR
line_search = none

petsc_options_iname=' -ksp_type            '   # -mat_view
petsc_options_value='  ksp_parrot_preonly  '   # ::ascii_matlab

dt = 1e7
num_steps=100

[./Quadrature] type=GRID order=TENTH [../]

[]

[Outputs]
file_base = AdvectionOut_${unifSteps}_${adapSteps}
exodus = true
csv=true
perf_graph = true
[]


[UserObjects]
[./soln]
type = SolveDiffusion
execute_on = 'initial'
block_id='1 6 4 7 2'
value_p ='1e-6 1e-6 1e-1 1e-6 1e-5'
boundary_D_bc = '6 7'
value_D_bc='4.0 1.0'
boundary_N_bc = ''
value_N_bc=''
aux_variable=P_aux
# output_file=matrix.e
[../]

[./MassAssembly]
type = AssembleMassMatrix
operator_userobject = storeOperatorsUO
block_id = '1 6 7 2 4'
value_p = ' 0.2 0.2 0.2 0.25 0.4'
#value_p = ' 1.0 1.0 1.0 1.0 1.0'
execute_on = 'initial'
constrain_matrix = true
dc_boundaries = '6'
dc_variables='CM'
value_D_bc='0.01'
[../]


[./storeOperatorsUO]
type = StoreOperators
[../]
[]
