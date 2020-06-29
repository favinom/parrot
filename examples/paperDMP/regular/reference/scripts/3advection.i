[GlobalParams]
conservative = false
[]

[Problem]
type = ParrotProblem
use_AFC = true
operator_userobject = storeOperatorsUO
[]

[Mesh]
file = ../refined_${mBe}_${mFe}.xdr
parallel_type = distributed
[]

[MeshModifiers]
[./fractureUserObject]
type = FractureUserObject
fn = 6
fx_string = '0.5,0.5,0.75,0.75,0.625,0.625'
fy_string = '0.5,0.5,0.75,0.75,0.625,0.625'
fa1_string = '0.0,90.0,0.0,90.0,0.0,90.0'
fd1_string = '1.0,1.0,0.5,0.5,0.25,0.25'
fd2_string = '1.0e-4,1.0e-4,1.0e-4,1.0e-4,1.0e-4,1.0e-4'
[../]
[]


[Variables]
[./CM] [../]
[]

[AuxVariables]
[./pressure] [../]
[./correction] [../]
[]

[Materials]
[./conductivity0]
type = FlowAndTransport k = 1.0 phi=1.0 block = 1 phiFrac=1.0 kFrac=1e4 pressure = pressure
fractureMeshModifier = fractureUserObject
[../] # conservative = true
[./conductivity1]
type = FlowAndTransport k = 1.0 phi=1.0 block = 2 phiFrac=1.0 kFrac=1e4 pressure = pressure
fractureMeshModifier = fractureUserObject
[../] # conservative = true
[]

[Kernels]
[upwind] type = Advection variable = CM [../]
[./time] type = PorosityTimeDerivative variable = CM lumping = true [../]
[]

[BCs]
[./u_injection_left] type = DirichletBC boundary = 'left' variable = CM value = 1.0 [../]
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

petsc_options_iname=' -ksp_type            ' # -mat_view
petsc_options_value='  ksp_parrot_preonly  ' # ::ascii_matlab

start_time = 0.0
end_time = 0.5
dt = 0.0025

[]

[Outputs]
file_base = AdvectionOut_${mBe}_${mFe}
[./exodus]
type = Exodus
sync_only = true
sync_times = '0.01 0.1 0.5'
[../]
csv = true
perf_graph = true
[]

[UserObjects]

[./soln]
type = SolveDiffusion3
execute_on = 'initial'
block_id='1 2'
value_p ='1 1 1e4'
boundary_D_bc='1'
value_D_bc='1.0'
boundary_N_bc='3 '
value_N_bc='-1.0 '
aux_variable=pressure
fractureMeshModifier = fractureUserObject
# output_file=DiffusionOut2R_${mBe}_${mFe}
solver_type = 1
conservative = false
[../]


[./storeOperatorsUO]
type = StoreOperators
[../]

[./MassAssembly]
type = AssembleMassMatrix
operator_userobject = storeOperatorsUO
block_id = '1 2'
value_p = '1.0 1.0'
execute_on = 'initial'
constrain_matrix = true
dc_variables='pressure'
dc_boundaries = '1'
value_D_bc='1.0'
#fractureMeshModifier = fractureUserObject
[../]

[]

