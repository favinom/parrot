[Problem]
type = ParrotProblem
[]

[Mesh]
file = refined_${resolution}_${unifSteps}_000${adaptSteps}_mesh.xdr
[]

[MeshModifiers]
[./aa2] type = BenchRegular2D [../]
[]

[Variables]
[./CM] [../]
[]
 
[AuxVariables]
[./P_aux] [../]
[]

[Materials]
[./conductivity0] type = FlowAndTransport k = 1.0 phi=1.0 block = 0 pressure = P_aux [../]
[./conductivity1] type = FlowAndTransport k = 1e4 phi=1.0 block = 1 pressure = P_aux [../]
[]
 
[AuxKernels]
[./en]
type = SolutionAux
solution = soln
variable = P_aux
scale_factor = 1.0
execute_on = 'initial'
[../]
[]

[Kernels]
active='time upwind'

[upwind]
type = AlgebraicDiffusion
variable = CM
[../]

[./time]
type = PorosityTimeDerivative
variable = CM
lumping = true
 dim = 2
[../]

[]

[BCs]
[./u_injection_left]
 type = DirichletBC boundary = 'left' variable = CM value = 1.0
[../]
[]

[Preconditioning]
[./SMP]
type = SMP
full = true
[../]
[]

 
 [Executioner]
 
 type = Transient
 solve_type= LINEAR
 line_search = none
 
 petsc_options_iname=' -ksp_type            ' # -mat_view
 petsc_options_value='  ksp_parrot_preonly  ' # ::ascii_matlab
 
 dt = 1.0e-3
 num_steps=500.0

 [./Quadrature]
 order=TENTH
 [../]

 
 []

 

[Outputs]
 file_base = AdvectionOut_${unifSteps}_${adaptSteps}
exodus = true
csv=true
perf_graph = true
[]


[UserObjects]
[./soln]
type = SolutionUserObject
mesh = DiffusionOut_${resolution}_${unifSteps}_${adaptSteps}.e
timestep = 2
system_variables = pressure
execute_on = 'initial'
[../]

# [./myuo]
# type = MyUO
# execute_on = 'linear'
# [../]
[]
