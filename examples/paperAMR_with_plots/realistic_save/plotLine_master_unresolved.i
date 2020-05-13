[Mesh]
 file = refinedMesh_${bepar}_0000_mesh.xdr
 parallel_type = DISTRIBUTED
[]

[Problem]
type = FEProblem
solve = false
kernel_coverage_check = false
[]

[AuxVariables]
[./P_aux]
[../]
[]

[UserObjects]
[./soln]
type = SolutionUserObject
# mesh = AdvectionOut_${bepar}_${amrpar}.e
mesh = DiffusionOut_${bepar}_${amrpar}.e
system_variables = pressure
execute_on = 'initial'
timestep = 'LATEST'
[../]
[]

[AuxKernels]
[./en]
type = SolutionAux
solution = soln
variable = P_aux
scale_factor = 1.0
execute_on = 'initial'
[../]
[]

[Executioner]
type = Transient
num_steps = 1
dt = 1

solve_type = 'PJFNK'

petsc_options_iname = '-pc_type -pc_hypre_type'
petsc_options_value = 'hypre boomeramg'
[]


[Transfers]

[./toSubHorizontal]
type = MultiAppProjectionTransfer
direction = to_multiapp
multi_app = sub_horizontal
source_variable = P_aux
variable = uh
[../]

[./toSubVertical]
type = MultiAppProjectionTransfer
direction = to_multiapp
multi_app = sub_vertical
source_variable = P_aux
variable = uv
[../]


[]


[MultiApps]

[./sub_horizontal]
type = TransientMultiApp
app_type = parrotApp
execute_on = timestep_end
input_files = plotLine_horizontal.i
cli_args = 'Outputs/out/file_base=horizontal_line_${bepar}_${amrpar}'
[../]

[./sub_vertical]
type = TransientMultiApp
app_type = parrotApp
execute_on = timestep_end
input_files = plotLine_vertical.i
cli_args = 'Outputs/out/file_base=vertical_line_${bepar}_${amrpar}'
[../]

[]
