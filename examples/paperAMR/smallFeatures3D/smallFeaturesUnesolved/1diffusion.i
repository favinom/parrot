[Mesh]
 file = refinedMesh_00${adapSteps}_mesh.xdr
  boundary_id = '11 22 23'
  boundary_name = 'inflow outflow1 outflow2'
# partitioner = linear
[]

[MeshModifiers]
[./fractureUserObject]
type = FractureUserObject
fn = 8
fx_string = '0.5,0.5,0.5,0.5,0.2,0.2,0.77,0.83'
fy_string = '1.125,0.175,1.6,1.6,2.05,2.05,2.05,2.05'
fz_string = '0.5,0.5,0.675,0.31,0.5,0.5,0.5,0.5'
fd1_string = '0.9,0.9,0.9,0.9,0.4,0.4,0.4,0.4'
fd2_string = '1.75,0.25,1.25,1.2472,0.30594,0.30594,0.3,0.3'
fd3_string = '0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01'
fa1_string = '0,0,0,0,78.6901,-78.6901,0,0'
fa2_string = '0,90,0,0,-90,-90,-90,-90'
fa3_string = '0,0,16.2602,-15.8192,90,-90,0,0'
[../]
[]


[Variables]
[./pressure] order=FIRST  family=LAGRANGE [../]
[]

[Kernels]
[./myDiffusion] type = PermeabilityDiffusion variable = pressure [../]
[]

[Materials]
[./conductivity1] 
type = FlowAndTransport
conservative = false
fractureMeshModifier =  fractureUserObject
phi = 0.0 phiFrac = 0.0
k = 1 kFrac = 1e4
[../]
[]

# observe that with the second BCs the stiffness matrix is SDP and we can use choleski factorization
[BCs]
[./inflowBC]   type = NeumannBC   variable = pressure value = 1.0  boundary = inflow  [../]
[./outflowBC1] type = DirichletBC variable = pressure value = 0.0  boundary = outflow1 [../]
[./outflowBC2] type = DirichletBC variable = pressure value = 0.0  boundary = outflow2 [../]
[]
 
[Preconditioning]
[./prec] type = SMP full = true ksp_norm = default [../]
[]
 
[Executioner]

 type=Steady
 solve_type= LINEAR
 line_search = none
 petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
 petsc_options_value='  preonly   lu       NONZERO               mumps '
 
# petsc_options_iname = '-pc_type -pc_hypre_type'
# petsc_options_value = 'hypre boomeramg'

[./Quadrature] order = NINTH type = GRID [../]
[]


[Outputs]
 file_base  = DiffusionOut_${adapSteps}
 exodus     = true
 perf_graph = true
[]
