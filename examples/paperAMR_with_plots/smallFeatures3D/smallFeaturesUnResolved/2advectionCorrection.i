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
[./conductivity1] type = FractureMaterial 
fractureMeshModifier =  fractureUserObject
matrixPorosity = 0.2 fracturePorosity = 0.2
matrixPermeability = 1 fracturePermeability = 1e4
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
#active='soln storeOperatorsUO'
[./soln]
type = SolveDiffusion
execute_on = 'initial'
block_id='0'
value_p ='1 1e4'
boundary_D_bc = '22 23'
value_D_bc='0.0 0.0'
boundary_N_bc = '11'
value_N_bc='-1.356070292717265' # '-1.3793251106'
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
block_id = '0 0'
value_p = ' 0.2 0.2'
execute_on = 'initial'
constrain_matrix = true
fractureMeshModifier = fractureUserObject
[../]
[./antidiffusiveFluxes]
 type = AntidiffusiveFluxes
operator_userobject = storeOperatorsUO
 execute_on = 'timestep_end'
 dc_boundaries = '11'
[../]
[]

[Postprocessors]
[./fluxBoundary]
  type = SideIntegralForFluxPostprocessor
  variable = P_aux
  boundary   = '11'
#  execute_on = 'initial'
[../]

[./Concentration0]
  type = ElementIntegralConcentrationPostprocessor
  variable = CM
  fractureRegionId = 0
  fractureMeshModifier =  fractureUserObject
#  execute_on = 'timestep_end'
[../]

[./volume0]
  type = ElementIntegralVolumePostprocessor
  fractureRegionId = 0
  fractureMeshModifier =  fractureUserObject
#  execute_on = 'initial'
[../]

[./Concentration1]
type = ElementIntegralConcentrationPostprocessor
variable = CM
fractureRegionId = 1
fractureMeshModifier =  fractureUserObject
#  execute_on = 'timestep_end'
[../]

[./volume1]
type = ElementIntegralVolumePostprocessor
fractureRegionId = 1
fractureMeshModifier =  fractureUserObject
#  execute_on = 'initial'
[../]
 
[./Concentration2]
 type = ElementIntegralConcentrationPostprocessor
 variable = CM
 fractureRegionId = 2
 fractureMeshModifier =  fractureUserObject
#  execute_on = 'timestep_end'
 [../]
 
[./volume2]
 type = ElementIntegralVolumePostprocessor
 fractureRegionId = 2
 fractureMeshModifier =  fractureUserObject
#  execute_on = 'initial'
[../]

[./Concentration3]
type = ElementIntegralConcentrationPostprocessor
variable = CM
fractureRegionId = 3
fractureMeshModifier =  fractureUserObject
#  execute_on = 'timestep_end'
[../]

[./volume3]
type = ElementIntegralVolumePostprocessor
fractureRegionId = 3
fractureMeshModifier =  fractureUserObject
#  execute_on = 'initial'
[../]

[./Concentration4]
type = ElementIntegralConcentrationPostprocessor
variable = CM
fractureRegionId = 4
fractureMeshModifier =  fractureUserObject
#  execute_on = 'timestep_end'
[../]

[./volume4]
type = ElementIntegralVolumePostprocessor
fractureRegionId = 4
fractureMeshModifier =  fractureUserObject
#  execute_on = 'initial'
[../]

[./Concentration5]
type = ElementIntegralConcentrationPostprocessor
variable = CM
fractureRegionId = 5
fractureMeshModifier =  fractureUserObject
#  execute_on = 'timestep_end'
[../]

[./volume5]
type = ElementIntegralVolumePostprocessor
fractureRegionId = 5
fractureMeshModifier =  fractureUserObject
#  execute_on = 'initial'
[../]

[./Concentration6]
type = ElementIntegralConcentrationPostprocessor
variable = CM
fractureRegionId = 6
fractureMeshModifier =  fractureUserObject
#  execute_on = 'timestep_end'
[../]

[./volume6]
type = ElementIntegralVolumePostprocessor
fractureRegionId = 6
fractureMeshModifier =  fractureUserObject
#  execute_on = 'initial'
[../]
 
[./Concentration7]
type = ElementIntegralConcentrationPostprocessor
variable = CM
fractureRegionId = 7
fractureMeshModifier =  fractureUserObject
#  execute_on = 'timestep_end'
[../]

[./volume7]
type = ElementIntegralVolumePostprocessor
fractureRegionId = 7
fractureMeshModifier =  fractureUserObject
#  execute_on = 'initial'
[../]
[]
