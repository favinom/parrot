[GlobalParams]
stabilize = false
[]

[Problem]
solve = false
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
[./pressure] order=FIRST  family=LAGRANGE [../]
[]
[AuxVariables]
[./flux_1][../]
[./flux_2][../]
[]

[Executioner]

type=Steady
solve_type= LINEAR
line_search = none
petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
petsc_options_value='  preonly   lu       NONZERO               mumps '

[./Quadrature] order = TENTH type=GRID [../]
[]



[Outputs]
file_base  = DiffusionOutS_${mResName}_${mRefLevName}_${mUmr}
exodus     = true
[]


[UserObjects]

[./soln]
type = SolveDiffusion3
execute_on = 'initial'
block_id='0 2 3 1'
value_p ='1e-6 1e-5 1e-6 1e-5 1e-1'
boundary_D_bc = '11 22'
value_D_bc='4.0 1.0'
boundary_N_bc = ''
value_N_bc=''
aux_variable='pressure'
fractureMeshModifier = fractureUserObject
solver_type = 1
conservative = false
[../]


[./storeOperatorsUO]
type = StoreOperators
[../]

[./MassAssembly]
type = AssembleMassMatrix
operator_userobject = storeOperatorsUO
block_id='0 2 3 1'
value_p = ' 0.2 0.25 0.2 0.25 0.4'
execute_on = 'initial'
constrain_matrix = true
dc_boundaries = '11'
dc_variables='pressure'
value_D_bc='0.01'
fractureMeshModifier = fractureUserObject
[../]


[./assF]
type = AssembleFlux
execute_on = 'timestep_end'
block_id='0 2 3 1'
value_p ='1e-6 1e-5 1e-6 1e-5 1e-1'
boundary_D_bc='11 22'
value_D_bc='4.0 1.0'
boundary_N_bc=''
value_N_bc=''
fractureMeshModifier = fractureUserObject
conservative = false
dc_boundaries = '11 22'
value_b=0.0
dc_variables='pressure'
sol_variable='pressure'
boundary_M_bc='new_side_set_1 new_side_set_2'
operator_userobject=storeOperatorsUO
[../]

[]
