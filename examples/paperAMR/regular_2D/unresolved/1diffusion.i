[GlobalParams]
conservative = true
[]

[Mesh]
 file = refinedMesh_00${adapSteps}_mesh.xdr
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

[Kernels]
[./myDiffusion] type = PermeabilityDiffusion variable = pressure [../]
[]

[Materials]
[./conductivity1] type = FlowAndTransport fractureMeshModifier =  fractureUserObject
phi = 1.0 phiFrac = 1.0
k = 1 kFrac = 1e4
[../]
[]

# observe that with the second BCs the stiffness matrix is SDP and we can use choleski factorization
[BCs]
[./inflowBC]  type = DirichletBC variable = pressure value = 1.0  boundary = right [../]
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
 
[./Quadrature] order = NINTH type = GRID [../]
[]


[Outputs]
 file_base  = DiffusionOut_${adapSteps}
 exodus     = true
 perf_graph = true
[]
