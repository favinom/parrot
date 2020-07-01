[GlobalParams]
conservative = false
stabilize = false
[]

[Problem]
type = ParrotProblem3
use_AFC = true
antidiffusive_fluxes=antidiffusiveFluxes
operator_userobject = storeOperatorsUO
solver_type=1
solve=false
[]

[Mesh]
file = Refined_${mRes}_0001_mesh.xdr
parallel_type=distributed
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

[./my]
type = FractureRefinement
fractureMeshModifier = fractureUserObject
refinements='${mRefLev} ${mUmr}'
doBoundaryRefinement = false
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
block = '0 3'
k = 1e-6 phi =0.2
pressure = P_aux
conservative = false
fractureMeshModifier =  fractureUserObject
kFrac=1e-1    phiFrac=0.4
[../]
[./conductivity2]
type = FlowAndTransport
block = '1 2'
k = 1e-5 phi =0.25
pressure = P_aux
conservative = false
fractureMeshModifier =  fractureUserObject
kFrac=1e-1    phiFrac=0.4
[../]
[]
[Kernels]
[upwind] type = Advection variable = CM [../]
[./time] type = PorosityTimeDerivative variable = CM lumping = true [../]
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
file_base = AdvectionOut_${mResName}_${mRefLevName}_${mUmr}
[./exodus]
type = Exodus
sync_only = true
sync_times = '1e8 5e8 1e9'
[../]
exodus=true
csv = true
perf_graph = true
[]

[UserObjects]
active='soln MassAssembly storeOperatorsUO antidiffusiveFluxes'
[./soln]
type = SolveDiffusion3
execute_on = 'initial'
block_id='0 2 3 1'
value_p ='1e-6 1e-5 1e-6 1e-5 1e-1'
boundary_D_bc = '11 22'
value_D_bc='4.0 1.0'
boundary_N_bc = ''
value_N_bc=''
aux_variable=P_aux
fractureMeshModifier = fractureUserObject
#output_file=matrix.e
conservative=false
solver_type=1
[../]

[./MassAssembly]
type = AssembleMassMatrix
operator_userobject = storeOperatorsUO
block_id='0 2 3 1'
value_p = ' 0.2 0.25 0.2 0.25 0.4'
execute_on = 'initial'
constrain_matrix = true
dc_boundaries = '11'
dc_variables='CM'
value_D_bc='0.01'
fractureMeshModifier = fractureUserObject
[../]



[./storeOperatorsUO]
type = StoreOperators
[../]

[./antidiffusiveFluxes]
type = AntidiffusiveFluxes3
operator_userobject = storeOperatorsUO
execute_on = 'timestep_end'
dc_boundaries = '11'
WriteCorrection=true
[../]
[]

