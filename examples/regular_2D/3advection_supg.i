[Problem]
type = ParrotProblem
[]

[Mesh]
file = refined_${resolution}_${unifSteps}_000${adapSteps}_mesh.xda
[]

[AuxVariables]
[./P_aux] [../]
[./vel_x]
family=MONOMIAL [../]
[./vel_y]
family=MONOMIAL [../]
[./vel_z]
family=MONOMIAL [../]
[./k0]
order=CONSTANT
family=MONOMIAL [../]
[./k1]
order=CONSTANT
family=MONOMIAL [../]
[./phi]
order=CONSTANT
family=MONOMIAL [../]
[]

[AuxKernels]
[./en]
type = SolutionAux
solution = soln
variable = P_aux
scale_factor = 1.0
execute_on = 'initial'
[../]
[vel_x]
type= CondactivityAux
variable=vel_x
comp_i=0
[../]
[vel_y]
type=CondactivityAux
variable=vel_y
comp_i=1
[../]
[./k0]
type = MaterialRealTensorValueAux
variable = k0
property = conductivityTensor [../]
[./k1]
type = MaterialRealTensorValueAux
variable = k1
property = conductivityTensor [../]
[]

[Variables]
[./CM]
[../]
[]

[Materials]
 active='conductivity1'
[./conductivity1]
 type =  HydraulicConductivity2D
 pressure=P_aux
 cond0=true
 K_matrix=1.0
 phi_m=1.0
 phi_f=1.0
 fn = 6
 fx_string = '0.5,0.5,0.75,0.75,0.625,0.625'
 fy_string = '0.5,0.5,0.75,0.75,0.625,0.625'
 fa1_string = '0.0,90.0,0.0,90.0,0.0,90.0'
 fd1_string = '1.0,1.0,0.5,0.5,0.25,0.25'
 fd2_string = '1e-4,1e-4,1e-4,1e-4,1e-4,1e-4'
 [../]
 []

[Kernels]
[./time]
 type = PorosityTimeDerivative
 variable = CM
 lumping = true
 dim = 2
 [../]
 
[./diff] type = Advection variable = CM int_by_parts = false [../]

[./advectionsupg]
 type = AdvectionSUPG
 coef=0.1
 use_h=true
 variable = CM
 p=P_aux
 [../]

 
 []

[Functions]
[./fun_n]
type = ParsedFunction
value = '1*(x<0.2500)*(y<0.2500)*(z<0.2500)'
[../]
[]

[BCs]
[./u_injection_left]
type = DirichletBC
boundary = 'left'
variable = CM
value = 1.0
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

petsc_options_iname=' -ksp_type            '
petsc_options_value='  ksp_parrot_preonly  '

dt = 1.0e-3
num_steps=100

[]

[Outputs]
file_base = AdvectionOutput_supg_01_${resolution}_${unifSteps}_000${adapSteps}
exodus = true
csv=true
print_perf_log = true
[]


[UserObjects]
[./soln]
type = SolutionUserObject
mesh = DiffusionOutput_${resolution}_${unifSteps}_000${adapSteps}.e
timestep = LATEST
system_variables = pressure
execute_on = 'initial'
[../]
[]





