[GlobalParams]
conservative = false
stabilize = false
[]

[Problem]
type = ParrotProblem3
use_AFC = true
solver_type = 1
antidiffusive_fluxes=antidiffusive_fluxes
operator_userobject = storeOperatorsUO
solve = false
[]

[Mesh]
type = GeneratedMesh
xmin= 0.0
xmax= 1.0
ymin= 0.0
ymax= 1.0
nx = ${mRes}
ny = ${mRes}
dim = 2
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

[./aaa]
type = FractureRefinement
fractureMeshModifier = fractureUserObject
refinements='${mRefLev} ${mUmr}'
# outputFileName = mesh_${mResName}_${mRefLevName}_${mUmr}.e
doBoundaryRefinement = false
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
[./conductivity1] type = FlowAndTransport fractureMeshModifier =  fractureUserObject
phi = 1.0 phiFrac = 1.0
k = 1 kFrac = 1e4
pressure = pressure
[../]
[]

[Kernels]
[upwind] type = Advection variable = CM [../]
[./time] type = PorosityTimeDerivative variable = CM lumping = true [../]
[]

[BCs]
[./u_injection_left]
 type = DirichletBC
 boundary = left
 variable = CM
 value='1.0' [../]
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

start_time = 0.0
end_time = 0.5
dt = 0.0025


[./Quadrature] order= NINTH type = GRID [../]

[]


[Outputs]
file_base = AdvectionOut_${mResName}_${mRefLevName}_${mUmr}
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
block_id='0'
value_p ='1.0 1e4'
boundary_D_bc='1'
value_D_bc='1.0'
boundary_N_bc='3 '
value_N_bc='-1.0 '
aux_variable=pressure
fractureMeshModifier = fractureUserObject
# output_file=DiffusionOut2_${mResName}_${mRefLevName}_${mUmr}.e
solver_type = 1
conservative = false
[../]

[./storeOperatorsUO]
type = StoreOperators
[../]

[./MassAssembly]
type = AssembleMassMatrix
operator_userobject = storeOperatorsUO
block_id = '0'
value_p = '1.0'
execute_on = 'initial'
constrain_matrix = true
dc_variables='pressure'
dc_boundaries = '3'
value_D_bc='1.0'
#fractureMeshModifier = fractureUserObject
[../]
 
[./antidiffusive_fluxes]
 type = AntidiffusiveFluxes3
 operator_userobject = storeOperatorsUO
 execute_on = 'timestep_end'
 dc_boundaries = '3'
 WriteCorrection=true
[../]
[]
