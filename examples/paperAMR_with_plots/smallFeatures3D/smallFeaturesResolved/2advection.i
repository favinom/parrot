[Problem]
type = ParrotProblem
use_AFC = true
operator_userobject = storeOperatorsUO
antidiffusive_fluxes = antidiffusiveFluxes
[]


[Mesh]
file = refinedMesh_${typeMesh}_${origLevel}_${Uref}_000${adapSteps}_mesh.xdr
boundary_id = '21 22'
boundary_name = 'inflow outflow'
[]

[MeshModifiers]
active=''
[./fractureUserObject]
type = FractureUserObject
fn = 8
fx_string = '0.5  ,0.5  ,0.77,0.83, 0.2,0.2  , 0.5,0.5'
fy_string = '1.125,0.175,2.05,2.05, 2.05,2.05 , 1.6,1.6'
fz_string = '0.5  ,0.5  ,0.5 ,0.5, 0.5,0.5 , 0.675,0.31'
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
[./correction] [../]
[]

[Materials]

[./porosity_1]
type = FlowAndTransport
conservative=true
block = '11 12 13'
k = 1.0
phi = 0.2
pressure=P_aux
[../]

[./porosity_2]
type = FlowAndTransport
conservative=false
block = '1 2 3 4 5 6 7 8'
k = 1.0e4
phi = 0.2
pressure=P_aux
[../]

#[./conductivity1] type = FractureMaterial fractureMeshModifier =  fractureUserObject
#matrixPorosity = 0.2 fracturePorosity = 0.2
#matrixPermeability = 1 fracturePermeability = 1e4
#pressure = P_aux
#[../]
[]

[Kernels]
active='time upwind'
[upwind] type = Advection variable = CM [../]
[./time] type = PorosityTimeDerivative variable = CM lumping = true [../]
[]

[BCs]
[./u_injection_left] type = DirichletBC boundary = 21 variable = CM value='1.0' [../]
[]

[Preconditioning]
[./SMP]
type = SMP
full = true
[../]
[]

[Executioner]

[./TimeIntegrator]
type = VoidIntegrator
[../]

type = Transient
solve_type= LINEAR
line_search = none

petsc_options_iname=' -ksp_type            '   # -mat_view
petsc_options_value='  ksp_parrot_preonly  '   # ::ascii_matlab

dt = 0.01
num_steps=100

[./Quadrature] type=GRID order=TENTH [../]

[]

[Outputs]
file_base = AdvectionOut_${typeMesh}_${origLevel}_${Uref}_${adapSteps}
exodus = true
csv=true
perf_graph = true
[]


[UserObjects]
#active='ComputeStatistics'

[./antidiffusiveFluxes]
type = AntidiffusiveFluxes
execute_on = 'timestep_end'
dc_boundaries = '6'
operator_userobject = storeOperatorsUO
WriteCorrection=true
[../]

[./soln]
type = SolveDiffusion
execute_on = 'initial'
block_id='1 2 3 4 5 6 7 8 11 12 13'
value_p ='1e4 1e4 1e4 1e4 1e4 1e4 1e4 1e4 1 1 1'
boundary_D_bc = '22'
value_D_bc='0.0 0.0'
boundary_N_bc = '21'
value_N_bc='-1.0'
aux_variable=P_aux
conservative=true
[../]



[./MassAssembly]
type = AssembleMassMatrix
operator_userobject = storeOperatorsUO
block_id='1 2 3 4 5 6 7 8 11 12 13'
value_p ='0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2'
execute_on = 'initial'
constrain_matrix = true
dc_boundaries = '21'
dc_variables='CM'
value_D_bc='1.0'
[../]

[./assembleVolumeVectors]
type=AssembleVolumeVectors
FractureRegions=true
NRegions=8
execute_on = 'initial'
block_id='1 2 3 4 5 6 7 8'
#fractureMeshModifier = fractureUserObject
[../]

[./storeOperatorsUO]
type = StoreOperators
[../]
[]


[Postprocessors]
[./c0] type = IntegralSolutionOverRegionFast region = 0 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./c1] type = IntegralSolutionOverRegionFast region = 1 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./c2] type = IntegralSolutionOverRegionFast region = 2 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./c3] type = IntegralSolutionOverRegionFast region = 3 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./c4] type = IntegralSolutionOverRegionFast region = 4 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./c5] type = IntegralSolutionOverRegionFast region = 5 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./c6] type = IntegralSolutionOverRegionFast region = 6 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./c7] type = IntegralSolutionOverRegionFast region = 7 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]

[./reg0] type = IntegralSolutionOverRegionFast region = 0 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg1] type = IntegralSolutionOverRegionFast region = 1 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg2] type = IntegralSolutionOverRegionFast region = 2 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg3] type = IntegralSolutionOverRegionFast region = 3 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg4] type = IntegralSolutionOverRegionFast region = 4 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg5] type = IntegralSolutionOverRegionFast region = 5 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg6] type = IntegralSolutionOverRegionFast region = 6 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg7] type = IntegralSolutionOverRegionFast region = 7 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]

[]
