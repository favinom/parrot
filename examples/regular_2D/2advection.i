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
[./CM] [../]
[]
 
[AuxVariables]
[./pressure] [../]
[]

[Materials]
[./conductivity0] type = FlowAndTransport k = 1.0 phi=1.0 block = 0 pressure = pressure [../]
[./conductivity1] type = FlowAndTransport k = 1e4 phi=1.0 block = 1 pressure = pressure [../]
[./hhhh] type = RegionMaterial regionMeshModifier = aa2 [../]
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

[./Quadrature] order=TENTH [../]
[]

[Outputs]
file_base = AdvectionOut_${resolution}_${unifSteps}_${adaptSteps}
exodus = true
csv=true
perf_graph = true
[]

[UserObjects]
[./myuo]
 type = SolveDiffusion
 execute_on = 'initial'
 block_id='0 1'
 value_p='1.0 1.0e4'
 boundary_D_bc='1'
 value_D_bc='1.0'
 boundary_N_bc='3 '
 value_N_bc='-1.0 '
 aux_variable='pressure'
 [../]
 []
