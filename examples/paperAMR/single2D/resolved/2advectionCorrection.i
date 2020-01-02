[Problem]
type = ParrotProblem
use_AFC = true
operator_userobject = storeOperatorsUO
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
[./correction] [../]
[]

[Materials]
[./conductivity1] type = FlowAndTransport k = 1e-6 phi=0.2  conservative = true block = '1 6 7' pressure = P_aux [../]
[./conductivity4] type = FlowAndTransport k = 1e-1 phi=0.4  conservative = true block = 4       pressure = P_aux [../]
[./conductivity2] type = FlowAndTransport k = 1e-5 phi=0.25 conservative = true block = 2       pressure = P_aux [../]
[]

[Kernels]
active='time upwind'
[upwind] type = Advection variable = CM [../]
[./time] type = PorosityTimeDerivative variable = CM lumping = true [../]
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

 petsc_options_iname=' -ksp_type            '   # -mat_view
 petsc_options_value='  ksp_parrot_preonly  '   # ::ascii_matlab

dt = 1e7
num_steps=100

[./Quadrature] order=TENTH [../]

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
conservative = true
# output_file=matrix.e
[../]
[./storeOperatorsUO]
type = StoreOperators
[../]
[./MassAssembly]
type = AssembleMassMatrix
operator_userobject = storeOperatorsUO 
block_id = '1 6 7 2 4'
value_p = ' 0.2 0.2 0.2 0.4 0.25'
execute_on = 'initial'
constrain_matrix = true
[../]
[./antidiffusiveFluxes]
 type = AntidiffusiveFluxes
 execute_on = 'timestep_end'
 dc_boundaries = '6 7'
operator_userobject = storeOperatorsUO
 [../]
[]
