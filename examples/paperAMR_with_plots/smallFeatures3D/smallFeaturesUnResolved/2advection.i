[GlobalParams]
conservative = false
[]

[Problem]
 type = ParrotProblem3
 use_AFC = true
 operator_userobject = storeOperatorsUO
 solver_type = 3
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
fx_string = '0.5,0.5,0.5,0.5,0.2,0.2,0.77,0.83'
fy_string = '1.125,0.175,1.6,1.6,2.05,2.05,2.05,2.05'
fz_string = '0.5,0.5,0.675,0.31,0.5,0.5,0.5,0.5'
fd1_string = '0.9,0.9,0.9,0.9,0.4,0.4,0.4,0.4'
fd2_string = '1.75,0.25,1.25,1.2472,0.30594,0.30594,0.3,0.3'
fd3_string = '0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01'
fa1_string = '0,0,0,0,78.6901,-78.6901,0,0'
fa2_string = '0,90,0,0,-90,-90,-90,-90'
fa3_string = '0,0,16.2602,-15.8192,90,-90,0,0'
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
type = FlowAndTransport fractureMeshModifier =  fractureUserObject
phi = 0.2 phiFrac = 0.2
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
num_steps=100

[./Quadrature] order = NINTH type = GRID [../]
[]

[Outputs]
 file_base = AdvectionOut_${adapSteps}
 exodus = true
 csv=true
 perf_graph = true
[]


[UserObjects]
[./soln]
type = SolveDiffusion2
execute_on = 'initial'
block_id='0'
value_p ='1 1e4'
boundary_D_bc = '22 23'
value_D_bc='0.0 0.0'
boundary_N_bc = '11'
value_N_bc='-1.0'
aux_variable=P_aux
fractureMeshModifier = fractureUserObject
#output_file=matrix.e
solver_type = 3
[../]

[./assembleVolumeVectors]
type=AssembleVolumeVectors
FractureRegions=false
NRegions=8
execute_on = 'initial'
#block_id='1 2 3 4 5 6 7 8'
FractureMeshModifier= fractureUserObject
[../]

[./storeOperatorsUO]
type = StoreOperators
[../]


[./MassAssembly]
type = AssembleMassMatrix
operator_userobject = storeOperatorsUO 
block_id = '0'
value_p = ' 0.2 0.2'
execute_on = 'initial'
constrain_matrix = true
fractureMeshModifier = fractureUserObject
dc_boundaries = '11'
dc_variables='CM'
value_D_bc='1.0'

[../]

[]


[Postprocessors]
#active=''
 [./int0] type = IntegralSolutionOverRegionFast region = 0 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
 [./int1] type = IntegralSolutionOverRegionFast region = 1 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
 [./int2] type = IntegralSolutionOverRegionFast region = 2 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
 [./int3] type = IntegralSolutionOverRegionFast region = 3 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
 [./int4] type = IntegralSolutionOverRegionFast region = 4 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
 [./int5] type = IntegralSolutionOverRegionFast region = 5 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
 [./int6] type = IntegralSolutionOverRegionFast region = 6 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
 [./int7] type = IntegralSolutionOverRegionFast region = 7 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]

 [./reg0] type = IntegralSolutionOverRegionFast region = 0 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
 [./reg1] type = IntegralSolutionOverRegionFast region = 1 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
 [./reg2] type = IntegralSolutionOverRegionFast region = 2 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
 [./reg3] type = IntegralSolutionOverRegionFast region = 3 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
 [./reg4] type = IntegralSolutionOverRegionFast region = 4 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
 [./reg5] type = IntegralSolutionOverRegionFast region = 5 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
 [./reg6] type = IntegralSolutionOverRegionFast region = 6 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
 [./reg7] type = IntegralSolutionOverRegionFast region = 7 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]


[]
