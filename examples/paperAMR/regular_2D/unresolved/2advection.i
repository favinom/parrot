[GlobalParams]
conservative = true
[]

[Problem]
type = ParrotProblem
use_AFC = true
#operator_userobject = storeOperatorsUO
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
#active='soln'
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
[../]
[./antidiffusiveFluxes]
 type = AntidiffusiveFluxes
operator_userobject = storeOperatorsUO
 execute_on = 'timestep_end'
 dc_boundaries = '1'
[../]
[]

[Postprocessors]
[./cm0]
type = ElementIntegralConcentrationPostprocessor
variable = CM fractureRegionId = 0 fractureMeshModifier =  aa2
[../]
[./cm1]
type = ElementIntegralConcentrationPostprocessor
variable = CM fractureRegionId = 1 fractureMeshModifier =  aa2
[../]
[./cm2]
type = ElementIntegralConcentrationPostprocessor
variable = CM fractureRegionId = 2 fractureMeshModifier =  aa2
[../]
[./cm3]
type = ElementIntegralConcentrationPostprocessor
variable = CM fractureRegionId = 4 fractureMeshModifier =  aa2
[../]
[./cm4]
type = ElementIntegralConcentrationPostprocessor
variable = CM fractureRegionId = 4 fractureMeshModifier =  aa2
[../]
[./cm5]
type = ElementIntegralConcentrationPostprocessor
variable = CM fractureRegionId = 5 fractureMeshModifier =  aa2
[../]
[./cm6]
type = ElementIntegralConcentrationPostprocessor
variable = CM fractureRegionId = 6 fractureMeshModifier =  aa2
[../]
[./cm7]
type = ElementIntegralConcentrationPostprocessor
variable = CM fractureRegionId = 7 fractureMeshModifier =  aa2
[../]
[./cm8]
type = ElementIntegralConcentrationPostprocessor
variable = CM fractureRegionId = 8 fractureMeshModifier =  aa2
[../]
[./cm9]
type = ElementIntegralConcentrationPostprocessor
variable = CM fractureRegionId = 9 fractureMeshModifier =  aa2
[../]

[./vol0]
type = ElementIntegralVolumePostprocessor
fractureRegionId = 0 fractureMeshModifier =  aa2
[../]
[./vol1]
type = ElementIntegralVolumePostprocessor
fractureRegionId = 1 fractureMeshModifier =  aa2
[../]
[./vol2]
type = ElementIntegralVolumePostprocessor
fractureRegionId = 2 fractureMeshModifier =  aa2
[../]
[./vol3]
type = ElementIntegralVolumePostprocessor
fractureRegionId = 3 fractureMeshModifier =  aa2
[../]
[./vol4]
type = ElementIntegralVolumePostprocessor
fractureRegionId = 4 fractureMeshModifier =  aa2
[../]
[./vol5]
type = ElementIntegralVolumePostprocessor
fractureRegionId = 5 fractureMeshModifier =  aa2
[../]
[./vol6]
type = ElementIntegralVolumePostprocessor
fractureRegionId = 6 fractureMeshModifier =  aa2
[../]
[./vol7]
type = ElementIntegralVolumePostprocessor
fractureRegionId = 7 fractureMeshModifier =  aa2
[../]
[./vol8]
type = ElementIntegralVolumePostprocessor
fractureRegionId = 8 fractureMeshModifier =  aa2
[../]
[./vol9]
type = ElementIntegralVolumePostprocessor
fractureRegionId = 9 fractureMeshModifier =  aa2
[../]

[]