[Problem]
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

[./subdomains_1]
type = SubdomainBoundingBox
top_right =   ' 0.501 1.001 0.0'
bottom_left = '-0.01 -0.01 0.0'
block_id = 1
block_name = block_1
[../]

[./subdomains_2]
type = SubdomainBoundingBox
top_right =   '1.01 1.01 0.0'
bottom_left = '0.499 -0.01 0.0'
block_id = 2
block_name = block_2
[../]


[./center_side_set_2]
type = SideSetsBetweenSubdomains
master_block = block_1
paired_block = block_2
new_boundary = 'new_side_set_2'
depends_on = 'subdomains_1 subdomains_2'
[../]

[./fractureUserObject]
type = FractureUserObject
fn = 6
fx_string = '0.5,0.5,0.75,0.75,0.625,0.625'
fy_string = '0.5,0.5,0.75,0.75,0.625,0.625'
fa1_string = '0.0,90.0,0.0,90.0,0.0,90.0'
fd1_string = '1.0,1.0,0.5,0.5,0.25,0.25'
fd2_string = '1.0e-4,1.0e-4,1.0e-4,1.0e-4,1.0e-4,1.0e-4'
[../]
[./my]
type = FractureRefinement
fractureMeshModifier = fractureUserObject
refinements='${mRefLev} ${mUmr}'
outputFileName = mesh_${mResName}_${mRefLevName}_${mUmr}.e
doBoundaryRefinement = true
[../]
[]


[Variables]
[./CM] order=FIRST  family=LAGRANGE [../]
[]

[AuxVariables]
[./P_aux] [../]
[./flux_1] order=FIRST  family=LAGRANGE [../]
[./flux_2] order=FIRST  family=LAGRANGE [../]
[]


[Preconditioning]
[./prec] type = SMP full = true ksp_norm = default [../]
[]

[Executioner]

type=Steady
solve_type= NEWTON
line_search = none
petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
petsc_options_value='  preonly   lu       NONZERO               mumps '

# petsc_options_iname = '-pc_type -pc_hypre_type'
# petsc_options_value = 'hypre boomeramg'

[./Quadrature] order = NINTH type = GRID [../]
[]

[UserObjects]
active=' '

[./soln]
type = SolveDiffusion3
execute_on = 'initial'
block_id='0 2'
value_p ='1 1 1e4'
boundary_D_bc='1'
value_D_bc='1.0'
boundary_N_bc='3 '
value_N_bc='-1.0 '
aux_variable=CM
fractureMeshModifier = fractureUserObject
output_file=DiffusionOut2_${mResName}_${mRefLevName}_${mUmr}.e
solver_type = 1
conservative = false
[../]

[./storeOperatorsUO]
type = StoreOperators
[../]

[./MassAssembly]
type = AssembleMassMatrix
operator_userobject = storeOperatorsUO
block_id = '0 2'
value_p = '1.0 1.0'
execute_on = 'initial'
constrain_matrix = true
dc_variables='CM'
dc_boundaries = '1'
value_D_bc='1.0'
#fractureMeshModifier = fractureUserObject
[../]


[./assF]
type = AssembleFlux
execute_on = 'timestep_end'
block_id='0 2'
value_p ='1.0 1.0 1e4'
boundary_D_bc='1'
#value_D_bc='1.0'
boundary_N_bc='3 '
value_N_bc='-1.0 '
fractureMeshModifier = fractureUserObject
conservative = false
value_b=0.0
operator_userobject = storeOperatorsUO
dc_variables='CM'
dc_boundaries = '1'
boundary_M_bc='new_side_set_2'
[../]

[]

[Outputs]
file_base = refinedMesh
perf_graph = true
[]
