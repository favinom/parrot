[Mesh]
file = refined_${resolution}_${unifSteps}_000${adaptSteps}_mesh.xdr
[]

[MeshModifiers]
[./aa2] type = BenchRegular2D [../]
[]

[Variables]
[./pressure] [../]
[]

[Materials]
[./conductivity0] type = FlowAndTransport k = 1.0 phi=1.0 block = 0 conservative = true [../]
[./conductivity1] type = FlowAndTransport k = 1e4 phi=1.0 block = 1 conservative = true [../]
[]

[Kernels]
[./time] type = PermeabilityDiffusion variable = pressure [../]
[]

[BCs]
[./right] type = DirichletBC boundary = 'right' variable = pressure value = 1.0 [../]
[./left] type = NeumannBC boundary = 'left' variable = pressure value = 1.0 [../]
[]

[Preconditioning]
[./SMP]
type = SMP
full = true
[../]
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
file_base = DiffusionOut_${resolution}_${unifSteps}_${adaptSteps}
exodus = true
perf_graph = true
[]
