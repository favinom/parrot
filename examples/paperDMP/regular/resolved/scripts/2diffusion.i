[GlobalParams]
stabilize = false
[]

[Problem]
solve = false
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

[./Quadrature] order = NINTH type = GRID [../]
[]



[Outputs]
file_base  = DiffusionOutR_${mBe}_${mFe}
exodus     = true
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


[./assF]
type = AssembleFlux3
execute_on = 'timestep_end'
block_id='1 2'
value_p ='1.0 1.0 1e4'
boundary_D_bc='1'
#value_D_bc='1.0'
boundary_N_bc='3 '
value_N_bc='-1.0 '
fractureMeshModifier = fractureUserObject
conservative = false
value_b=0.0
operator_userobject = storeOperatorsUO
dc_variables='pressure'
dc_boundaries = '1'
sol_variable='pressure'
boundary_M_bc='new_side_set_2'
[../]

[]
