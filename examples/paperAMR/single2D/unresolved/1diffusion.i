[Mesh]
 file = refinedMesh_00${adapSteps}_mesh.xdr
  block_id = '1 6 4 7 2'
  boundary_id = '11 22'
  boundary_name = 'inflow outflow'
# partitioner = linear
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
[./conductivity2] type = FractureMaterial fractureMeshModifier =  fractureUserObject block = 2
matrixPorosity = 0.0 fracturePorosity = 0.0
matrixPermeability = 1e-5 fracturePermeability = 1e-1
[../]
[]

# observe that with the second BCs the stiffness matrix is SDP and we can use choleski factorization
[BCs]
[./inflowBC]  type = DirichletBC variable = pressure value = 4.0  boundary = inflow  [../]
[./outflowBC] type = DirichletBC variable = pressure value = 1.0  boundary = outflow [../]
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
