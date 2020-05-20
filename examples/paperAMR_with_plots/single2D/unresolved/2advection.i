[Problem]
type = ParrotProblem
use_AFC = true
operator_userobject = storeOperatorsUO
antidiffusive_fluxes=antidiffusive_fluxes
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
[./correction] [../]
[]

[Materials]
[./conductivity1]
type = FlowAndTransport
block = 0
k = 1 phi =2.0e5
pressure = P_aux
conservative = false
fractureMeshModifier =  fractureUserObject
kFrac=1e5    phiFrac=4e5
[../]
[./conductivity2]
type = FlowAndTransport
block = 2
k = 10 phi =2.5e5
pressure = P_aux
conservative = false
fractureMeshModifier =  fractureUserObject
kFrac=1e5    phiFrac=4e5
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
file_base = AdvectionOut1000_${adapSteps}
exodus = true
csv=true
perf_graph = true
[]


[UserObjects]

[./soln]
type = SolveDiffusion
execute_on = 'initial'
block_id='0 2'
value_p ='1 10 1e5'
boundary_D_bc = '11 22'
value_D_bc='4.0 1.0'
boundary_N_bc = ''
value_N_bc=''
aux_variable=P_aux
fractureMeshModifier = fractureUserObject
#output_file=matrix.e
conservative=false
[../]

[./MassAssembly]
type = AssembleMassMatrix
operator_userobject = storeOperatorsUO
block_id = '0   2'
value_p = ' 2e5 2.5e5 4e4'
execute_on = 'initial'
constrain_matrix = true
dc_boundaries = '11'
dc_variables='CM'
value_D_bc='0.01'
fractureMeshModifier = fractureUserObject
[../]
 
[./antidiffusive_fluxes]
 type = AntidiffusiveFluxes
 operator_userobject = storeOperatorsUO
 execute_on = 'timestep_end'
 dc_boundaries = '11 22'
 WriteCorrection=false
[../]

[./storeOperatorsUO]
type = StoreOperators
[../]
[]










