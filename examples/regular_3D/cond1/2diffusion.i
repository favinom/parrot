[Problem]
type = ParrotProblem
use_AFC = false
[]

[Mesh]
file = ../refined_${resolution}_${unifSteps}_000${adaptSteps}_mesh.xdr
[]

[MeshModifiers]
[./aa2] type = RegionsRegular3D [../]
[]

[Variables]
[./pressure] [../]
[]

[Materials]
[./conductivity0] type = FlowAndTransport k = 1.0 phi=0.1 block = 0 [../]
[./conductivity1] type = FlowAndTransport k = 1.0e-4 phi=0.01 block = 1 [../]
[./conductivity2] type = FlowAndTransport k = 0.1 phi=0.1 block = 2 [../]
[]

[Kernels]
[./time] type = PermeabilityDiffusion variable = pressure [../]
[]

[BCs]
[./dirBC]  type = DirichletBC variable = pressure value = 1.0  boundary = '14' [../]
[./fluxBC] type = NeumannBC   variable = pressure value = 1.0  boundary = '13' [../]
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
