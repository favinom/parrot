[GlobalParams]
conservative = false
[]

[Problem]
type = ParrotProblem
use_AFC = true
antidiffusive_fluxes=antidiffusive
operator_userobject = storeOperatorsUO
[]

[Mesh]
 file = refinedMesh_${resolution}_00${adapSteps}_mesh.xdr
[]

[MeshModifiers]
[./fractureUserObject]
 type = FractureUserObject
 fn = 6
 fx_string = '0.5,0.5,0.75,0.75,0.625,0.625'
 fy_string = '0.5,0.5,0.75,0.75,0.625,0.625'
 fa1_string = '0.0,90.0,0.0,90.0,0.0,90.0'
 fd1_string = '1.0,1.0,0.5,0.5,0.25,0.25'
 fd2_string = '1.0e-4,1.0e-4,1.0e-4,1.0e-4,1.0e-4,1.0e-4'
[../]
[./aa2] type = BenchRegular2D [../]
[]

[Variables]
[./CM] [../]
[]
 
[AuxVariables]
[./P_aux] [../]
[./correction] [../]
[]

[Materials]
[./conductivity1] type = FlowAndTransport fractureMeshModifier =  fractureUserObject
phi = 1.0 phiFrac = 1.0
k = 1 kFrac = 1e4
pressure = P_aux
[../]
[]

[Kernels]
active='time upwind'
[upwind] type = Advection variable = CM [../]
[./time] type = PorosityTimeDerivative variable = CM lumping = true [../]
[]

[BCs]
[./u_injection_left] type = DirichletBC boundary = left variable = CM value='1.0' [../]
[]

[Preconditioning]
[./SMP]
type = SMP
full = true
[../]
[]

[Executioner]

[./TimeIntegrator] type = VoidIntegrator [../]

type = Transient
solve_type= LINEAR
line_search = none

 petsc_options_iname=' -ksp_type            '   # -mat_view
 petsc_options_value='  ksp_parrot_preonly  '   # ::ascii_matlab

 dt = 0.0025
 num_steps=100.0

[./Quadrature] order= NINTH type = GRID [../]

[]

[Outputs]
 file_base = AdvectionOut_${resolution}_${adapSteps}
 exodus = true
 csv=true
 perf_graph = true
[]


[UserObjects]
 
[./antidiffusive]
 type = AntidiffusiveFluxes
 operator_userobject = storeOperatorsUO
 dc_boundaries = '1'
 execute_on = 'timestep_end'
 WriteCorrection=true
[../]
 
[./soln]
type = SolveDiffusion
execute_on = 'initial'
block_id='0'
value_p ='1 1e4'
boundary_D_bc='1'
value_D_bc='1.0'
boundary_N_bc='3 '
value_N_bc='-1.0 '
aux_variable=P_aux
fractureMeshModifier = fractureUserObject
#output_file=matrix.e
[../]
 
[./storeOperatorsUO]
type = StoreOperators
[../]
 
[./MassAssembly]
type = AssembleMassMatrix
operator_userobject = storeOperatorsUO 
block_id = '0'
value_p = ' 1.0 1.0'
execute_on = 'initial'
constrain_matrix = true
fractureMeshModifier = fractureUserObject
dc_boundaries = '3'
dc_variables='CM'
value_D_bc='1.0'
[../]

 
[./assembleVolumeVectors]
type=AssembleVolumeVectors
fractureMeshModifier = aa2
execute_on = 'initial'
[../]
[]

[Postprocessors]
active=''
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
