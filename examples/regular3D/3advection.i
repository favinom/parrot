[GlobalParams]
conservative = true
[]

[Problem]
type = ParrotProblem
use_AFC = true
#operator_userobject = storeOperatorsUO
[]

[Mesh]
file = refined_${resolution}_${unifSteps}_000${adaptSteps}_mesh.xdr
[]

#[MeshModifiers]
#[./aa2] type = BenchRegular2D [../]
#[]

[Variables]
[./CM] [../]
[]
 
[AuxVariables]
[./pressure] [../]
[./correction] [../]
[]

[Materials]
[./conductivity0] type = FlowAndTransport k = 1.0   phi=0.1   block = 0 pressure = pressure [../]
[./conductivity1] type = FlowAndTransport k = 0.1   phi=0.1   block = 1 pressure = pressure [../]
[./conductivity2] type = FlowAndTransport k = ${kf} phi=${pf} block = 2 pressure = pressure [../]
[]

[Kernels]
active='time upwind'
[upwind] type = AlgebraicDiffusion variable = CM [../]
[./time] type = PorosityTimeDerivative variable = CM lumping = true
[../]
[]

[BCs]
[./u_injection_left] type = DirichletBC boundary = 13 variable = CM value = 1.0 [../]
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
 
 dt = 0.0025
 num_steps=100.0

[./Quadrature] order=TENTH [../]
[]

[Outputs]
file_base = AdvectionOut_${resolution}_${unifSteps}_${adaptSteps}
exodus = true
csv=true
perf_graph = true
[]

[UserObjects]
#active = 'soln'
[./soln]
 type = SolveDiffusion
 execute_on = 'initial'
 block_id='0   1   2'
 value_p= '1.0 0.1 ${kf}'
 boundary_D_bc='10'
 value_D_bc='1.0'
 boundary_N_bc='13'
 value_N_bc='-1.0'
 aux_variable='pressure'
 output_file=matrix.e
#conservative = true
 [../]
[./storeOperatorsUO]
type = StoreOperators
[../]
[./MassAssembly]
type = AssembleMassMatrix
operator_userobject = storeOperatorsUO 
block_id = '0 1'
value_p = ' 1.0 1.0'
execute_on = 'initial'
constrain_matrix = true
[../]
[./antidiffusiveFluxes]
 type = AntidiffusiveFluxes
 execute_on = 'timestep_end'
 dc_boundaries = '13'
operator_userobject = storeOperatorsUO
 [../]
 []

[Postprocessors]
active=''
 [./int0] type = ElementIntegralVariableOverRegion variable = CM region = 0 [../]
 [./int1] type = ElementIntegralVariableOverRegion variable = CM region = 1 [../]
 [./int2] type = ElementIntegralVariableOverRegion variable = CM region = 2 [../]
 [./int3] type = ElementIntegralVariableOverRegion variable = CM region = 3 [../]
 [./int4] type = ElementIntegralVariableOverRegion variable = CM region = 4 [../]
 [./int5] type = ElementIntegralVariableOverRegion variable = CM region = 5 [../]
 [./int6] type = ElementIntegralVariableOverRegion variable = CM region = 6 [../]
 [./int7] type = ElementIntegralVariableOverRegion variable = CM region = 7 [../]
 [./int8] type = ElementIntegralVariableOverRegion variable = CM region = 8 [../]
 [./int9] type = ElementIntegralVariableOverRegion variable = CM region = 9 [../]
 
 [./reg0] type = ElementIntegralVariableOverRegion region = 0 [../]
 [./reg1] type = ElementIntegralVariableOverRegion region = 1 [../]
 [./reg2] type = ElementIntegralVariableOverRegion region = 2 [../]
 [./reg3] type = ElementIntegralVariableOverRegion region = 3 [../]
 [./reg4] type = ElementIntegralVariableOverRegion region = 4 [../]
 [./reg5] type = ElementIntegralVariableOverRegion region = 5 [../]
 [./reg6] type = ElementIntegralVariableOverRegion region = 6 [../]
 [./reg7] type = ElementIntegralVariableOverRegion region = 7 [../]
 [./reg8] type = ElementIntegralVariableOverRegion region = 8 [../]
 [./reg9] type = ElementIntegralVariableOverRegion region = 9 [../]
 []

