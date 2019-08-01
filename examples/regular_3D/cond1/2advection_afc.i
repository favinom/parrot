[Problem]
 type = ParrotProblem
 use_AFC = true
 []
 
[Mesh]
 file = ../refined_${resolution}_${unifSteps}_000${adaptSteps}_mesh.xdr
 []
 
[MeshModifiers]
[./aa2] type = RegionsRegular3D [../]
 []
 
[Variables]
[./CM] [../]
 []
 
[AuxVariables]
[./pressure] [../]
 []
 
[Materials]
[./conductivity0] type = FlowAndTransport k = 1.0 phi=0.1 block = 0 pressure = pressure [../]
[./conductivity1] type = FlowAndTransport k = 1.0e-4 phi=0.01 block = 1 pressure = pressure [../]
[./conductivity2] type = FlowAndTransport k = 0.1 phi=0.1 block = 2 pressure = pressure [../]
[./hhhh] type = RegionMaterial regionMeshModifier = aa2 [../]
[]

[Kernels]
 active='time upwind'
 
[upwind]
 type = Advection
 variable = CM
 [../]
 
[./time]
 type = PorosityTimeDerivative
 variable = CM
 lumping = true
 dim = 3
 [../]
 
 []
 
[BCs]
[./u_injection_left]
 type = DirichletBC boundary = '13' variable = CM value = 1.0
 [../]
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
[./myuo]
 type = SolveDiffusion
 execute_on = 'initial'
 block_id='0 1 2'
 value_p='1.0 1.0e4 0.1'
 boundary_D_bc='14'
 value_D_bc='1.0'
 boundary_N_bc='13 '
 value_N_bc='-1.0 '
 aux_variable='pressure'
 [../]
 []
 
[Postprocessors]
 
[./intC00]
 type = ElementIntegralVariableOverRegion
 variable = CM
 region = 0
 execute_on = 'timestep_end'
 [../]
 
[./intD00]
 type = ElementIntegralVariableOverRegion
 region = 0
 execute_on = 'timestep_end'
 [../]
 
[./intC01]
 type = ElementIntegralVariableOverRegion
 variable = CM
 region = 1 execute_on = 'timestep_end'
 [../]
[./intD01]
 type = ElementIntegralVariableOverRegion
 region = 1 execute_on = 'timestep_end'
 [../]
 
[./intC02]
 type = ElementIntegralVariableOverRegion
 variable = CM
 region = 2 execute_on = 'timestep_end'
 [../]
 
[./intD02]
 type = ElementIntegralVariableOverRegion
 region = 2 execute_on = 'timestep_end'
 [../]
 
[./intC03]
 type = ElementIntegralVariableOverRegion
 variable = CM
 region = 3 execute_on = 'timestep_end'
 [../]
 
[./intD03]
 type = ElementIntegralVariableOverRegion
 region = 3 execute_on = 'timestep_end'
 [../]
 
[./intC04]
 type = ElementIntegralVariableOverRegion
 variable = CM
 region = 4 execute_on = 'timestep_end'
 [../]
[./intD04]
 type = ElementIntegralVariableOverRegion
 region = 4 execute_on = 'timestep_end'
 [../]
 
[./intC05]
 type = ElementIntegralVariableOverRegion
 variable = CM
 region = 5 execute_on = 'timestep_end'
 [../]
[./intD05]
 type = ElementIntegralVariableOverRegion
 region = 5 execute_on = 'timestep_end'
 [../]
 
[./intC06]
 type = ElementIntegralVariableOverRegion
 variable = CM
 region = 6 execute_on = 'timestep_end'
 [../]
[./intD06]
 type = ElementIntegralVariableOverRegion
 region = 6 execute_on = 'timestep_end'
 [../]
 
[./intC07]
 type = ElementIntegralVariableOverRegion
 variable = CM
 region = 7 execute_on = 'timestep_end'
 [../]
[./intD07]
 type = ElementIntegralVariableOverRegion
 region = 7 execute_on = 'timestep_end'
 [../]
 
[./intC08]
 type = ElementIntegralVariableOverRegion
 variable = CM
 region = 8 execute_on = 'timestep_end'
 [../]
[./intD08]
 type = ElementIntegralVariableOverRegion
 region = 8 execute_on = 'timestep_end'
 [../]
 
[./intC09]
 type = ElementIntegralVariableOverRegion
 variable = CM
 region = 9 execute_on = 'timestep_end'
 [../]
[./intD09]
 type = ElementIntegralVariableOverRegion
 region = 9 execute_on = 'timestep_end'
 [../]
 
[./intC10]
 type = ElementIntegralVariableOverRegion
 variable = CM
 region = 10 execute_on = 'timestep_end'
 [../]
[./intD10]
 type = ElementIntegralVariableOverRegion
 region = 10 execute_on = 'timestep_end'
 [../]
 
[./intC11]
 type = ElementIntegralVariableOverRegion
 variable = CM
 region = 11 execute_on = 'timestep_end'
 [../]
[./intD11]
 type = ElementIntegralVariableOverRegion
 region = 11 execute_on = 'timestep_end'
 [../]
 
[./intC12]
 type = ElementIntegralVariableOverRegion
 variable = CM
 region = 12 execute_on = 'timestep_end'
 [../]
[./intD12]
 type = ElementIntegralVariableOverRegion
 region = 12 execute_on = 'timestep_end'
 [../]
 
[./intC13]
 type = ElementIntegralVariableOverRegion
 variable = CM
 region = 13 execute_on = 'timestep_end'
 [../]
[./intD13]
 type = ElementIntegralVariableOverRegion
 region = 13 execute_on = 'timestep_end'
 [../]
 
[./intC14]
 type = ElementIntegralVariableOverRegion
 variable = CM
 region = 14 execute_on = 'timestep_end'
 [../]
[./intD14]
 type = ElementIntegralVariableOverRegion
 region = 14 execute_on = 'timestep_end'
 [../]
 
[./intC15]
 type = ElementIntegralVariableOverRegion
 variable = CM
 region = 15 execute_on = 'timestep_end'
 [../]
[./intD15]
 type = ElementIntegralVariableOverRegion
 region = 15 execute_on = 'timestep_end'
 [../]
 
[./intC16]
 type = ElementIntegralVariableOverRegion
 variable = CM
 region = 16 execute_on = 'timestep_end'
 [../]
[./intD16]
 type = ElementIntegralVariableOverRegion
 region = 16 execute_on = 'timestep_end'
 [../]
 
[./intC17]
 type = ElementIntegralVariableOverRegion
 variable = CM
 region = 17 execute_on = 'timestep_end'
 [../]
[./intD17]
 type = ElementIntegralVariableOverRegion
 region = 17 execute_on = 'timestep_end'
 [../]
 
[./intC18]
 type = ElementIntegralVariableOverRegion
 variable = CM
 region = 18
 execute_on = 'timestep_end'
 [../]
 
[./intD18]
 type = ElementIntegralVariableOverRegion
 region = 18 execute_on = 'timestep_end'
 [../]
 
[./intC19]
 type = ElementIntegralVariableOverRegion
 variable = CM
 region = 19 execute_on = 'timestep_end'
 [../]
[./intD19]
 type = ElementIntegralVariableOverRegion
 region = 19 execute_on = 'timestep_end'
 [../]
 
[./intC20]
 type = ElementIntegralVariableOverRegion
 variable = CM
 region = 20 execute_on = 'timestep_end'
 [../]
[./intD20]
 type = ElementIntegralVariableOverRegion
 region = 20 execute_on = 'timestep_end'
 [../]
 
[./intC21]
 type = ElementIntegralVariableOverRegion
 variable = CM
 region = 21 execute_on = 'timestep_end'
 [../]
[./intD21]
 type = ElementIntegralVariableOverRegion
 region = 21 execute_on = 'timestep_end'
 [../]
 
 []
