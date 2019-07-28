[Problem]
type = ParrotProblem
use_AFC = false
[]

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
[./conductivity0] type = FlowAndTransport k = 1.0 phi=1.0 block = 0 [../]
[./conductivity1] type = FlowAndTransport k = 1e4 phi=1.0 block = 1 [../]
[]

[Kernels]
[./time] type = PermeabilityDiffusion variable = pressure [../]
[]

[BCs]
[./right]
 type = DirichletBC boundary = 'right' variable = pressure value = 1.0
[../]
[./left]
 type = NeumannBC boundary = 'left' variable = pressure value = 1.0
 [../]

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
 
 petsc_options_iname=' -ksp_type            ' # -mat_view
 petsc_options_value='  ksp_parrot_preonly  ' # ::ascii_matlab
 
[./Quadrature] order=TENTH [../]
[]

[Outputs]
file_base = DiffusionOut_${resolution}_${unifSteps}_${adaptSteps}
exodus = true
perf_graph = true
[]
