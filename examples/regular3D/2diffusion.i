#[Problem]
#type = ParrotProblem
#[]

[Mesh]
file = refined_${resolution}_${unifSteps}_000${adaptSteps}_mesh.xdr
[]

[MeshModifiers]
[./aa1]
type = FractureUserObject
fn = 9
fx_string = '0.5,0.5,0.5,0.749975,0.75,0.749975,0.625,0.625,0.625'
fy_string = '0.5,0.5,0.5,0.749975,0.749975,0.75,0.625,0.625,0.625'
fz_string = '0.5,0.5,0.5,0.75,0.749975,0.749975, 0.625,0.625,0.625'
fa1_string = '0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0'
fa2_string = '0.0,90.0,0.0,0.0,90.0,0.0,0.0,90.0,0.0'
fa3_string = '0.0,0.0,90.0,0.0,0.0,90.0,0.0,0.0,90.0'
fd1_string = '1.0,1.0,1.0,0.50005,0.50005,0.50005,0.2501,0.2501,0.2501'
fd2_string = '1.0,1.0,1.0,0.50005,0.50005,0.50005,0.2501,0.2501,0.2501'
fd3_string = '0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001'
[../]
[]

[Variables]
[./pressure] [../]
[]

[Materials]
[./conductivity0] type = FlowAndTransport k = 1.0 phi=1.0 conservative = true block = 0 [../]
[./conductivity1] type = FlowAndTransport k = 0.1 phi=1.0 conservative = true block = 1 [../]
[./conductivity2] type = FlowAndTransport k = 1e4 phi=1.0 conservative = true block = 2 [../]
[]

[Kernels]
[diffusion] type =  PermeabilityDiffusion variable = pressure [../]
[]

[BCs]
[./dirBC] type = DirichletBC variable = pressure value = 1 boundary = '10' [../]
[./neuBC] type = NeumannBC   variable = pressure value = 1 boundary = '13' [../]
[]

[Preconditioning]
[./SMP]
type = SMP
full = true
[../]
[]

 
 [Executioner]
 
 type = Steady
 solve_type= LINEAR
 line_search = none
 
 petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
 petsc_options_value='  preonly   lu       NONZERO               mumps '
 
# petsc_options_iname = '-pc_type -pc_hypre_type'
# petsc_options_value = 'hypre boomeramg'
 
[./Quadrature] order=TENTH [../]
[]

[Outputs]
file_base = DiffusionOut_${resolution}_${unifSteps}_${adaptSteps}
exodus = true
perf_graph = true
[]

[UserObjects]
[./soln]
 type = ComputeStatistics
 execute_on = 'initial'
 fractureMeshModifier = aa1
 [../]
[]