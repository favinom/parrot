[Problem]
type = ParrotProblem
use_AFC = true
operator_userobject = storeOperatorsUO
[]

[Mesh]
file = refinedMesh_00${adapSteps}_mesh.xdr
boundary_id = '11 22 23'
boundary_name = 'inflow outflow1 outflow2'
[]

[MeshModifiers]
[./fractureUserObject]
type = FractureUserObject
fn = 8
fx_string =  '0.5  ,0.5  ,0.77,0.83, 0.2,0.2  , 0.5,0.5'
fy_string =  '1.125,0.175,2.05,2.05, 2.05,2.05 , 1.6,1.6'
fz_string =  '0.5  ,0.5  ,0.5 ,0.5, 0.5,0.5 , 0.675,0.31'
fa1_string = '0,90,90,90,78.6901,-78.6901,0,0'
fa2_string = '0, 0, 0, 0,0,0,0,0'
fa3_string = '0,90,90,90,90,90,16.2602,-15.8192'
fd1_string = '0.9,0.25,0.3,0.3,0.3059,0.3059,0.9,0.9'
fd2_string = '1.75,0.9,0.4,0.4,0.4,0.4,1.25,1.2472'
fd3_string = '0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01'
[../]
[]

[Variables]
[./CM] [../]
[]

[AuxVariables]
[./P_aux] [../]
[]

[Materials]
[./conductivity1]
type = FractureMaterial
fractureMeshModifier =  fractureUserObject
matrixPorosity = 0.2
fracturePorosity = 0.2
matrixPermeability = 1
fracturePermeability = 1e4
pressure = P_aux
conservative = false
[../]
[]

[Kernels]
active='time upwind'
[upwind] type = Advection variable = CM [../]
[./time] type = PorosityTimeDerivative variable = CM lumping = true [../]
[]

[BCs]
[./u_injection_left] type = DirichletBC boundary = 11 variable = CM value='1' [../]
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

dt = 0.01
num_steps=1

[./Quadrature]
order = TENTH type = GRID
[../]
[]

[Outputs]
file_base = AdvectionOut_${adapSteps}
exodus = true
csv=true
perf_graph = true
[]


[UserObjects]
[./soln]
type = SolveDiffusion
execute_on = 'initial'
block_id='0'
value_p ='1 1e4'
boundary_D_bc = '22 23'
value_D_bc='0.0 0.0'
boundary_N_bc = '11'
value_N_bc='-1.3793251106'
#'-1.3492491344'
#'-1.3379067369'
aux_variable=P_aux
fractureMeshModifier = fractureUserObject
conservative=false
[../]
 
[./MassAssembly]
type = AssembleMassMatrix
operator_userobject = storeOperatorsUO
fractureMeshModifier = fractureUserObject
block_id = '0'
value_p = ' 0.2 0.2'
execute_on = 'initial'
constrain_matrix = true
dc_boundaries = '11'
dc_variables='CM'
value_D_bc='1.0'
[../]
 
[./storeOperatorsUO]
 type = StoreOperators
[../]
 
[./assembleVolumeVectors]
type=AssembleVolumeVectors
FractureRegions=true
NRegions=8
execute_on = 'initial'
#block_id='1 2 3 4 5 6 7 8'
fractureMeshModifier = fractureUserObject
[../]
[]


[Postprocessors]
[./fluxBoundary]
type = SideIntegralForFluxPostprocessor
variable = P_aux
boundary   = '11'
#  execute_on = 'initial'
[../]




[./volume0]
type = ElementIntegralVolumePostprocessor
fractureRegionId = 0
fractureMeshModifier =  fractureUserObject
#  execute_on = 'initial'
[../]



[./volume1]
type = ElementIntegralVolumePostprocessor
fractureRegionId = 1
fractureMeshModifier =  fractureUserObject
#  execute_on = 'initial'
[../]



[./volume2]
type = ElementIntegralVolumePostprocessor
fractureRegionId = 2
fractureMeshModifier =  fractureUserObject
#  execute_on = 'initial'
[../]



[./volume3]
type = ElementIntegralVolumePostprocessor
fractureRegionId = 3
fractureMeshModifier =  fractureUserObject
#  execute_on = 'initial'
[../]


[./volume4]
type = ElementIntegralVolumePostprocessor
fractureRegionId = 4
fractureMeshModifier =  fractureUserObject
#  execute_on = 'initial'
[../]



[./volume5]
type = ElementIntegralVolumePostprocessor
fractureRegionId = 5
fractureMeshModifier =  fractureUserObject
#  execute_on = 'initial'
[../]



[./volume6]
type = ElementIntegralVolumePostprocessor
fractureRegionId = 6
fractureMeshModifier =  fractureUserObject
#  execute_on = 'initial'
[../]


[./volume7]
type = ElementIntegralVolumePostprocessor
fractureRegionId = 7
fractureMeshModifier =  fractureUserObject
#  execute_on = 'initial'
[../]

 [./reg0] type = IntegralSolutionOverRegionFast region = 0 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
 [./reg1] type = IntegralSolutionOverRegionFast region = 1 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
 [./reg2] type = IntegralSolutionOverRegionFast region = 2 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
 [./reg3] type = IntegralSolutionOverRegionFast region = 3 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
 [./reg4] type = IntegralSolutionOverRegionFast region = 4 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
 [./reg5] type = IntegralSolutionOverRegionFast region = 5 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
 [./reg6] type = IntegralSolutionOverRegionFast region = 6 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
 [./reg7] type = IntegralSolutionOverRegionFast region = 7 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
 
[]

