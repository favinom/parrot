[Mesh]
  file = refinedMesh_00${adapSteps}_mesh.xdr
  boundary_id = '11 22'
  boundary_name = 'inflow outflow'
  parallel_type=distributed
#  uniform_refine=3
# partitioner = linear
[]

[MeshModifiers]
[./fractureUserObject]
type = FractureUserObject
fn = 1
fx_string = '0.0'
fy_string = '0.0'
fa1_string = '-30.963756532073525'
fd1_string = '2000.0'
fd2_string = '0.01'
[../]
[]


[Variables]
[./CM] order=FIRST  family=LAGRANGE [../]
[]
[AuxVariables]
[./P_aux] [../]
[./flux_1][../]
[./flux_2][../]
[./correction] [../]
[]
 
[Kernels]
[./myDiffusion] type = PermeabilityDiffusion variable = CM [../]
[]

[Materials]
[./conductivity1]
 type = FlowAndTransport
 block = '0 3'
 k = 1e-6 phi =0.2
 pressure = P_aux
 conservative = false
 fractureMeshModifier =  fractureUserObject
 kFrac=1e-1    phiFrac=0.4
 [../]
[./conductivity2]
 type = FlowAndTransport
 block = '1 2'
 k = 1e-5 phi =0.25
 pressure = P_aux
 conservative = false
 fractureMeshModifier =  fractureUserObject
 kFrac=1e-1    phiFrac=0.4
 [../]
 []

# observe that with the second BCs the stiffness matrix is SDP and we can use choleski factorization
[BCs]
[./inflowBC]  type = DirichletBC variable = CM value = 4.0  boundary = inflow  [../]
[./outflowBC] type = DirichletBC variable = CM value = 1.0  boundary = outflow [../]
[]
 
[Preconditioning]
[./prec] type = SMP full = true ksp_norm = default [../]
[]
 
[Executioner]

 type=Steady
 solve_type= LINEAR
 line_search = none
 petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
 petsc_options_value='  preonly   lu       NONZERO               mumps '
 
[./Quadrature] order = TENTH type=GRID [../]
[]


[Outputs]
 file_base  = DiffusionOut_${adapSteps}
 exodus     = true
 perf_graph = true
[]

[UserObjects]
 active='assF storeOperatorsUO MassAssembly'
[./soln]
 type = SolutionUserObject
 mesh = AdvectionOut_80_08.e
 timestep = LATEST
 system_variables = P_aux
 execute_on = 'initial'
 [../]
 
[./assF]
 type = AssembleFlux
 execute_on = 'timestep_end'
 block_id='0 2 3 1'
 value_p ='1e-6 1e-5 1e-6 1e-5 1e-1'
 boundary_D_bc='11 22'
 value_D_bc='4.0 1.0'
 boundary_N_bc=''
 value_N_bc=''
 fractureMeshModifier = fractureUserObject
 conservative = false
 dc_boundaries = '11 22'
 value_b=0.0
 dc_variables='CM'
 boundary_M_bc='new_side_set_1 new_side_set_2'
 operator_userobject=storeOperatorsUO
 [../]
 
[./MassAssembly]
 type = AssembleMassMatrix
 operator_userobject = storeOperatorsUO
 block_id = '1 6 7 2 4'
 value_p = ' 0.2 0.2 0.2 0.25 0.4'
#value_p = ' 1.0 1.0 1.0 1.0 1.0'
 execute_on = 'initial'
 constrain_matrix = true
 dc_boundaries = '6'
 dc_variables='CM'
 value_D_bc='0.01'
 [../]
 
 [./antidiffusive_fluxes]
 type = AntidiffusiveFluxes
 operator_userobject = storeOperatorsUO
 execute_on = 'timestep_end'
 dc_boundaries = '11 22'
 WriteCorrection=false
 [../]
[./storeOperatorsUO]
 type = StoreOperators
 [../]
 []
