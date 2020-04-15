[Mesh]
file = mesh_${typeMesh}_${origLevel}.e
boundary_id = '21 22 23'
boundary_name = 'inflow outflow1 outflow2'
uniform_refine = ${Uref}
# partitioner = linear
[]

[MeshModifiers]
 active=''
[./fractureUserObject]
type = FractureUserObject
fn = 0
fx_string = '0.5  ,0.5  ,0.77,0.83, 0.2,0.2  , 0.5,0.5'
fy_string = '1.125,0.175,2.05,2.05, 2.05,2.05 , 1.6,1.6'
fz_string = '0.5  ,0.5  ,0.5 ,0.5, 0.5,0.5 , 0.675,0.31'
fa1_string = '0,90,90,90,78.6901,-78.6901,0,0'
fa2_string = '0, 0, 0, 0,0,0,0,0'
fa3_string = '0,90,90,90,90,90,16.2602,-15.8192'
fd1_string = '0.9,0.25,0.3,0.3,0.3059,0.3059,0.9,0.9'
fd2_string = '1.75,0.9,0.4,0.4,0.4,0.4,1.25,1.2472'
fd3_string = '0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01'
[../]
[]


[Variables]
[./pressure] order=FIRST  family=LAGRANGE [../]
[]

[Kernels]
[./myDiffusion] type = PermeabilityDiffusion variable = pressure [../]
[]


[Materials]

active='porosity_1 porosity_2'

[./conductivity1]
type = FlowAndTransport fractureMeshModifier =  fractureUserObject
phi = 0.2 phiFrac = 0.2
k = 1 kFrac = 1e4
pressure = pressure
conservative=false
[../]


[./porosity_1]
type = FlowAndTransport
conservative=false
block = '11 12 13 2 3 4 5 6 7 8'
k = 1.0
phi = 0.2
pressure=pressure
[../]

[./porosity_2]
type = FlowAndTransport
conservative=false
block = '1'
k = 1.0e2
phi = 0.2
pressure=pressure
[../]

[./porosity_3]
type = FlowAndTransport
conservative=false
block = '1'
k = 1.0
phi = 0.2
pressure=pressure
[../]
[]
# observe that with the second BCs the stiffness matrix is SDP and we can use choleski factorization
[BCs]
[./inflowBC]   type = NeumannBC   variable = pressure value = 1.0  boundary = 21 [../]
[./outflowBC1] type = DirichletBC variable = pressure value = 0.0  boundary = 22 [../]
[./outflowBC2] type = DirichletBC variable = pressure value = 0.0  boundary = 23 [../]
[]

[Preconditioning]
[./prec] type = SMP full = true ksp_norm = default [../]
[]

[Executioner]

type=Steady
solve_type= LINEAR
line_search = none
#petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
#petsc_options_value='  preonly   lu       NONZERO               mumps '

petsc_options_iname = '-pc_type -pc_hypre_type'
petsc_options_value = 'hypre boomeramg'

[./Quadrature] type=GRID order = TENTH [../]
[]


[Outputs]
file_base  = DiffusionOut_${typeMesh}_${origLevel}_${Uref}.e
exodus     = true
perf_graph = true
[]
