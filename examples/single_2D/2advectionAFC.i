[Problem]
type = ParrotProblem
use_AFC = true
change_sol=true
operator_userobject_problem = operator_userobject_problem
dc_boundaries = "6"
[]

[Mesh]
 file = refinedMesh_${unifSteps}_000${adapSteps}_mesh.xdr
 block_id = '1 6 4 7 2'
 boundary_id = '6 7'
 boundary_name = 'inflow outflow'
 uniform_refine = 0
#partitioner = linear
 []

[Variables]
[./CM] [../]
[]
 
[AuxVariables]
[./P_aux] [../]
[./correction] [../]
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
type = Advection
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
num_steps=100

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
[./operator_userobject_problem]
type = StoreOperators
#execute_on = 'initial'
[../]
[./MassAssembly]
type = AssembleMassMatrix
operator_userobject = operator_userobject_problem 
block_id = '1 6 7 2 4'
value_p = ' 0.2 0.2 0.2 0.4 0.25'
execute_on = 'initial'
[../]
[./PrintAssembly]
 type = PrintMatrix
 execute_on = 'timestep_end'
 dc_boundaries = "6"
 [../]
[]
