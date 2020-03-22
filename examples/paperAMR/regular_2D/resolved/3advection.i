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

[MeshModifiers]
[./aa2] type = BenchRegular2D [../]
[]

[Variables]
[./CM] [../]
[]
 
[AuxVariables]
[./pressure] [../]
[./correction] [../]
[]

[Materials]
[./conductivity0] type = FlowAndTransport k = 1.0 phi=1.0 block = 0 pressure = pressure [../] # conservative = true
[./conductivity1] type = FlowAndTransport k = 1e4 phi=1.0 block = 1 pressure = pressure [../] # conservative = true
[./hhhh] type = RegionMaterial regionMeshModifier = aa2 [../]
[]

[Kernels]
active='time upwind'
[upwind] type = AlgebraicDiffusion variable = CM [../]
[./time] type = PorosityTimeDerivative variable = CM lumping = true
[../]
[]

[BCs]
[./u_injection_left] type = DirichletBC boundary = 'left' variable = CM value = 1.0 [../]
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
 block_id='0 1'
 value_p='1.0 1.0e4'
 boundary_D_bc='1'
 value_D_bc='1.0'
 boundary_N_bc='3 '
 value_N_bc='-1.0 '
 aux_variable='pressure'
# output_file=matrix.e
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
 dc_boundaries = '1./'
operator_userobject = storeOperatorsUO
 [../]
[./assembleVolumeVectors]
type=AssembleVolumeVectors
fractureMeshModifier = aa2
execute_on = 'initial'
[../]
[]

[Postprocessors]

 [./int0] type = IntegralSolutionOverRegionFast region = 0 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
 [./int1] type = IntegralSolutionOverRegionFast region = 1 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
 [./int2] type = IntegralSolutionOverRegionFast region = 2 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
 [./int3] type = IntegralSolutionOverRegionFast region = 3 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
 [./int4] type = IntegralSolutionOverRegionFast region = 4 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
 [./int5] type = IntegralSolutionOverRegionFast region = 5 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
 [./int6] type = IntegralSolutionOverRegionFast region = 6 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
 [./int7] type = IntegralSolutionOverRegionFast region = 7 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
 [./int8] type = IntegralSolutionOverRegionFast region = 8 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
 [./int9] type = IntegralSolutionOverRegionFast region = 9 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]

 [./reg0] type = IntegralSolutionOverRegionFast region = 0 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
 [./reg1] type = IntegralSolutionOverRegionFast region = 1 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
 [./reg2] type = IntegralSolutionOverRegionFast region = 2 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
 [./reg3] type = IntegralSolutionOverRegionFast region = 3 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
 [./reg4] type = IntegralSolutionOverRegionFast region = 4 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
 [./reg5] type = IntegralSolutionOverRegionFast region = 5 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
 [./reg6] type = IntegralSolutionOverRegionFast region = 6 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
 [./reg7] type = IntegralSolutionOverRegionFast region = 7 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
 [./reg8] type = IntegralSolutionOverRegionFast region = 8 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
 [./reg9] type = IntegralSolutionOverRegionFast region = 9 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]

 []

