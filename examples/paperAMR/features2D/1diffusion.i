[Mesh]
 file = refinedMesh_00${adapSteps}_mesh.xdr
# partitioner = linear
[]

[MeshModifiers]
[./fractureUserObject]
type = FractureUserObject
fn = 7
fx_string = '0.0,50,-60,-50,-100.0,-125.0,-150.0'
fy_string = '100.0,-100,-100,0.0,25.0,50.0,-75.0'
fa1_string = '90.0,-63.4349488229,-116.5650511771,0.0,90.0,0.0,90.0'
fd1_string = '200.0,223.60679775,223.60679775,100.0,50.0,50.0,250.0'
fd2_string = '0.5,0.5,0.5,0.5,0.5,0.5,0.5'
[../]
[]

[Variables]
[./pressure] order=FIRST  family=LAGRANGE [../]
[]

[Kernels]
[./myDiffusion] type = PermeabilityDiffusion variable = pressure [../]
[]

[Materials]
[./conductivity1] type = FractureMaterial fractureMeshModifier =  fractureUserObject block = 0
matrixPorosity = 0.0 fracturePorosity = 0.0
matrixPermeability = 1e-6 fracturePermeability = 1e-1
[../]
[]

# observe that with the second BCs the stiffness matrix is SDP and we can use choleski factorization
[BCs]
[./inflowBC]  type = DirichletBC variable = pressure value = 4.0  boundary = top  [../]
[./outflowBC] type = DirichletBC variable = pressure value = 1.0  boundary = bottom [../]
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

[./Quadrature] order = TENTH [../]
[]


[Outputs]
 file_base  = DiffusionOut_${adapSteps}
 exodus     = true
 perf_graph = true
[]
