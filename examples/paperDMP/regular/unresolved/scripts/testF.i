[Problem]
type = FEProblem
solve = true
[]

[Mesh]
 file = 1refinedMesh_0001_mesh.xdr
 #mesh_0${mRes}_0${mRefLev}_0.xdr
[]
 
 
[Variables]
[./CM] [../]
[]
 
[AuxVariables]
    [./P_aux] [../]
    [./flux_1][../]
    [./flux_2][../]
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
 
[./left_block]
    type = SubdomainBoundingBox
    block_id = 1
    block_name = left_block
    bottom_left = '0.0 0.0 0.0'
    top_right =   '0.49995 1.0 0.0'
[../]
[./right_block]
    type = SubdomainBoundingBox
    block_id = 2
    block_name = right_block
    bottom_left = '0.4999500001 0.0 0.0'
    top_right = '1 1 0.0'
[../]
[./center_side_set]
    type = SideSetsBetweenSubdomains
    master_block = left_block
    paired_block = right_block
    new_boundary = 'new_side_set'
    depends_on = 'left_block right_block'
[../]

[]


[AuxKernels]
 active=''
[./en]
    type = SolutionAux
    solution = soln
    variable = P_aux
    scale_factor = 1.0
    execute_on = 'initial'
[../]
[]



[Materials]
[./conductivity1]
    type = FlowAndTransport
    fractureMeshModifier =  fractureUserObject
    phi = 1.0 phiFrac = 1.0
    k = 1.0 kFrac = 1.e4
    pressure = P_aux
    conservative=false
[../]
[]


[Kernels]
[./myDiffusion]
    type = PermeabilityDiffusion
    variable = CM
[../]
[]


[BCs]
[./inflowBC]
    type = DirichletBC
    variable = CM
    value = 1.0
    boundary = 1
[../]
[./ouflowBC]
    type = NeumannBC
    variable = CM
    value = 1.0
    boundary = 3
[../]
[]

[Preconditioning]
[./SMP]
    type = SMP
    full = true
[../]
[]

[Executioner]
    type=Steady
    solve_type= LINEAR
    line_search = none
    petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
    petsc_options_value='  preonly   lu       NONZERO               mumps '

  [./Quadrature] order = NINTH type = GRID [../]
[]

[Outputs]
    file_base = AdvectionOut_1280_7
    exodus = true
    csv=true
    print_perf_log = true
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
    type = AssembleFluxFracture
    execute_on = 'timestep_end'
    block_id='1 2'
    value_p ='1 1 1e4'
    boundary_D_bc='1'
    value_D_bc='1.0'
    boundary_N_bc='3 '
    boundary_M_bc='4 '
    value_N_bc='-1.0 '
    sol_variable=P_aux
    fractureMeshModifier = fractureUserObject
    conservative = false
    dc_boundaries = '1'
    matrix_1 = 1.0
    matrix_2 = 1.0
    dc_variables='CM'
    operator_userobject=storeOperatorsUO
[../]
[./storeOperatorsUO]
    type = StoreOperators
[../]
[./MassAssembly]
    type = AssembleMassMatrix
    operator_userobject = storeOperatorsUO
    block_id = '1 2'
    value_p = ' 1.0 1.0'
    execute_on = 'initial'
    constrain_matrix = true
    fractureMeshModifier = fractureUserObject
    value_D_bc='1.0'
    dc_variables='CM'
    dc_boundaries='1'
[../]
[]


