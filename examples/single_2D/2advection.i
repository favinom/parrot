[Problem]
type = ParrotProblem
use_AFC = false
[]

[Mesh]
 file = refinedMesh_${unifSteps}_000${adapSteps}_mesh.xdr
 block_id = '1 6 4 7 2'
 boundary_id = '6 7'
 boundary_name = 'inflow outflow'
 []

[Variables]
[./CM] [../]
[]
 
[AuxVariables]
[./P_aux] [../]
[]

[Materials]
[./conductivity1] type = FlowAndTransport k = 1e-6 phi=0.2 block = '1 6 7' pressure = P_aux [../]
[./conductivity4] type = FlowAndTransport k = 1e-1 phi=0.4 block = 4       pressure = P_aux [../]
[./conductivity2] type = FlowAndTransport k = 1e-5 phi=0.25 block = 2      pressure = P_aux [../]
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
[./u_injection_left] type = DirichletBC boundary = inflow variable = CM value='0.01' [../]
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

petsc_options_iname=' -ksp_type             '   # -mat_view
petsc_options_value='  ksp_parrot_preonly    '   # ::ascii_matlab

dt = 1e7
num_steps = 2.0 #100

[./Quadrature]
order=TENTH
[../]

[]

[Outputs]
 file_base = AdvectionOut_${unifSteps}_${adapSteps}
exodus = true
csv=true
perf_graph = true
[]


[UserObjects]
[./soln]
type = SolutionUserObject
mesh = DiffusionOut_${unifSteps}_${adapSteps}.e
timestep = LATEST
system_variables = pressure
execute_on = 'initial'
[../]

# [./myuo]
# type = MyUO
# execute_on = 'timestep_end'
# [../]
 []
