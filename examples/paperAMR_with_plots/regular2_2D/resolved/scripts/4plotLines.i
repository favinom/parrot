[Mesh]
 file = '../refined_${mBe}_${mFe}.xdr'
 parallel_type = DISTRIBUTED
[]

[Problem]
type = FEProblem
solve = false
kernel_coverage_check = false
[]

[AuxVariables]
[./P_aux][../]
[./CM1][../]
[./CM2][../]
[./CM3][../]
[]

[UserObjects]
[./pressureUO]
type = SolutionUserObject
mesh = '../AdvectionOut_${mBe}_${mFe}.e'
system_variables = pressure
execute_on = 'initial'
timestep = 'LATEST'
[../]
[./CM1UO]
type = SolutionUserObject
mesh = '../AdvectionOut_${mBe}_${mFe}.e'
system_variables = CM
execute_on = 'initial'
timestep = 1
[../]
[./CM2UO]
type = SolutionUserObject
mesh = '../AdvectionOut_${mBe}_${mFe}.e'
system_variables = CM
execute_on = 'initial'
timestep = 2
[../]
[./CM3UO]
type = SolutionUserObject
mesh = '../AdvectionOut_${mBe}_${mFe}.e'
system_variables = CM
execute_on = 'initial'
timestep = 3
[../]

[]

[AuxKernels]
[./pressureAK]
type = SolutionAux
solution = pressureUO
variable = P_aux
scale_factor = 1.0
execute_on = 'initial'
[../]
[./CM1AK]
type = SolutionAux
solution = CM1UO
variable = CM1
scale_factor = 1.0
execute_on = 'initial'
[../]
[./CM2AK]
type = SolutionAux
solution = CM2UO
variable = CM2
scale_factor = 1.0
execute_on = 'initial'
[../]
[./CM3AK]
type = SolutionAux
solution = CM3UO
variable = CM3
scale_factor = 1.0
execute_on = 'initial'
[../]

[]

[Executioner]
type = Transient
num_steps = 1
dt = 1
[]


[Transfers]
[./toPressureLineA]
type = MultiAppProjectionTransfer direction = to_multiapp multi_app = lineA source_variable = 'P_aux' variable = 'p' [../]
[./toCM1LineA]
type = MultiAppProjectionTransfer direction = to_multiapp multi_app = lineA source_variable = 'CM1' variable = 'cm1' [../]
[./toCM2LineA]
type = MultiAppProjectionTransfer direction = to_multiapp multi_app = lineA source_variable = 'CM2' variable = 'cm2' [../]
[./toCM3LineA]
type = MultiAppProjectionTransfer direction = to_multiapp multi_app = lineA source_variable = 'CM3' variable = 'cm3' [../]

[./toPressureLineB]
type = MultiAppProjectionTransfer direction = to_multiapp multi_app = lineB source_variable = 'P_aux' variable = 'p' [../]
[./toCM1LineB]
type = MultiAppProjectionTransfer direction = to_multiapp multi_app = lineB source_variable = 'CM1' variable = 'cm1' [../]
[./toCM2LineB]
type = MultiAppProjectionTransfer direction = to_multiapp multi_app = lineB source_variable = 'CM2' variable = 'cm2' [../]
[./toCM3LineB]
type = MultiAppProjectionTransfer direction = to_multiapp multi_app = lineB source_variable = 'CM3' variable = 'cm3' [../]

[./toPressureLineC]
type = MultiAppProjectionTransfer direction = to_multiapp multi_app = lineC source_variable = 'P_aux' variable = 'p' [../]
[./toCM1LineC]
type = MultiAppProjectionTransfer direction = to_multiapp multi_app = lineC source_variable = 'CM1' variable = 'cm1' [../]
[./toCM2LineC]
type = MultiAppProjectionTransfer direction = to_multiapp multi_app = lineC source_variable = 'CM2' variable = 'cm2' [../]
[./toCM3LineC]
type = MultiAppProjectionTransfer direction = to_multiapp multi_app = lineC source_variable = 'CM3' variable = 'cm3' [../]

[]


[MultiApps]

[./lineA]
type = TransientMultiApp
app_type = parrotApp
execute_on = timestep_end
input_files = plotLine_A.i
cli_args = 'Outputs/out/file_base=line_A_${mBe}_${mFe}'
[../]
[./lineB]
type = TransientMultiApp
app_type = parrotApp
execute_on = timestep_end
input_files = plotLine_B.i
cli_args = 'Outputs/out/file_base=line_B_${mBe}_${mFe}'
[../]
[./lineC]
type = TransientMultiApp
app_type = parrotApp
execute_on = timestep_end
input_files = plotLine_C.i
cli_args = 'Outputs/out/file_base=line_C_${mBe}_${mFe}'
[../]

[]

#[Outputs]
#exodus = true
#[]
