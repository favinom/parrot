[Mesh]
 file = refined_${resolution}_${unifSteps}_000${adaptSteps}_mesh.xdr
[]

[Variables]
[./pressure] order=FIRST  family=LAGRANGE [../]
[]

[Kernels]
[./myDiffusion] type = PermeabilityDiffusion variable = pressure [../]
[]

[Materials]
[./conductivity0] type = FlowAndTransport k = 1.0 phi=1.0 block = 0 [../]
[./conductivity1] type = FlowAndTransport k = 1e4 phi=1.0 block = 1 [../]
[]

[BCs]
[./dirBC]  type = DirichletBC variable = pressure value = 1  boundary = 'right'  [../]
[./fluxBC] type = NeumannBC variable = pressure value = 1  boundary = 'left' [../]
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

[./Quadrature]
 order = TENTH
[../]

[]


[Outputs]
 file_base  = DiffusionOut_${resolution}_${unifSteps}_${adaptSteps}
 exodus     = true
 perf_graph = true
[]

