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
[./pressure] order=FIRST  family=LAGRANGE [../]
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
[./soln]
type = SolveDiffusion3
execute_on = 'initial'
block_id='0'
value_p ='1 1e4'
boundary_D_bc='1'
value_D_bc='1.0'
boundary_N_bc='3 '
value_N_bc='-1.0 '
aux_variable=pressure
fractureMeshModifier = fractureUserObject
output_file=DiffusionOut2_${mResName}_${mRefLevName}_${mUmr}.e
solver_type = 1
conservative = false
[../]
[]

