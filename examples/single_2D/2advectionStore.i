[Problem]
type = ParrotProblem
[]

[Mesh]
 file = refinedMesh_000${adapSteps}_mesh.xdr
 block_id = '1 6 4 7 2'
 boundary_id = '6 7'
 boundary_name = 'inflow outflow'
 uniform_refine = 0
#second_order = true
 []

[Variables]
[./CM] [../]
[]
 
[AuxVariables]
[./P_aux] [../]
[]

[Materials]
[./conductivity1] type = FlowAndTransport k = 1e-6 phi=0.2 block = '1 6 7' [../]
[./conductivity4] type = FlowAndTransport k = 1e-1 phi=0.4 block = 4       [../]
[./conductivity2] type = FlowAndTransport k = 1e-5 phi=0.25 block = 2       [../]
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

[Materials]

[./epsInt2] type = GenericConstantMaterial block = 4 prop_names = epsInt prop_values = 0.01  [../]
[./epsInt4] type = GenericConstantMaterial block = 1 prop_names = epsInt prop_values = 0.0   [../]
[./epsInt5] type = GenericConstantMaterial block = 6 prop_names = epsInt prop_values = 0.0   [../]
[./epsInt6] type = GenericConstantMaterial block = 7 prop_names = epsInt prop_values = 0.0   [../]
[./epsInt7] type = GenericConstantMaterial block = 2 prop_names = epsInt prop_values = 0.0   [../]

[./dummy2] type = GenericConstantMaterial block = 4 prop_names = dummy prop_values = 0.0  [../]
[./dummy4] type = GenericConstantMaterial block = 1 prop_names = dummy prop_values = 0.0  [../]
[./dummy5] type = GenericConstantMaterial block = 6 prop_names = dummy prop_values = 0.0  [../]
[./dummy6] type = GenericConstantMaterial block = 7 prop_names = dummy prop_values = 0.0  [../]
[./dummy7] type = GenericConstantMaterial block = 2 prop_names = dummy prop_values = 1.0 [../]

[]


[Kernels]
active='time upwind'

[upwind]
type = AlgebraicDiffusion
#upwinding_type=full
variable = CM
#int_by_parts=false
[../]

[./time]
type = PorosityTimeDerivative
variable = CM
lumping = true
[../]

[]

[BCs]
[./u_injection_left] type = DirichletBC boundary = inflow variable = CM value='0.01' [../]
#[./u_injection_right] type = OutflowBC boundary = 'outflow 3' variable = CM [../]
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
order=SIXTH
[../]

[]

[Outputs]
 file_base = AdvectionOut_${adapSteps}
exodus = true
csv=true
print_perf_log = true
[]


[UserObjects]
[./soln]
type = SolutionUserObject
mesh = DiffusionOut_${adapSteps}.e
timestep = 2
system_variables = pressure
execute_on = 'initial'
[../]
[]

[Postprocessors]

[./flux_left]
type = SideFlux
variable = CM
boundary = outflow
# PER FAVORE CONTROLLARE IL COEF
coef = 0.0
execute_on = 'timestep_end'
[../]

[./int3]
type = ElementIntegral_phi_c_MatProp
variable = CM
mat_prop = dummy
execute_on = 'timestep_end'
[../]

[./intFrac]
type = ElementIntegral_phi_c_MatProp
variable = CM
mat_prop = epsInt
execute_on = 'timestep_end'
[../]

[]





