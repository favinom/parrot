[Mesh]
 file = refinedMesh_${unifSteps}_000${adapSteps}_mesh.xdr
  block_id = '1 6 4 7 2'
  boundary_id = '6 7'
  boundary_name = 'inflow outflow'
# partitioner = linear
[]

[Variables]
[./pressure] order=FIRST  family=LAGRANGE [../]
[]

[Kernels]
[./myDiffusion] type = PermeabilityDiffusion variable = pressure [../]
[]

[Materials]
[./conductivity1] type = FlowAndTransport k = 1e-6 conservative = true phi=0.0 block = '1 6 7' [../]
[./conductivity4] type = FlowAndTransport k = 1e-1 conservative = true phi=0.0 block = 4       [../]
[./conductivity2] type = FlowAndTransport k = 1e-5 conservative = true phi=0.0 block = 2       [../]
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
 file_base  = DiffusionOut_${unifSteps}_${adapSteps}
 exodus     = true
 perf_graph = true
[]
