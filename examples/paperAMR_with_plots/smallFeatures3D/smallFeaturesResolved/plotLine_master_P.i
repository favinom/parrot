[Mesh]
file = refinedMesh_${typeMesh}_${origLevel}_${Uref}_000${adapSteps}_mesh.xdr
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
#active='soln'
[./soln]
type = SolutionUserObject
mesh = AdvectionOut_${typeMesh}_${origLevel}_${Uref}_${adapSteps}.e
system_variables = P_aux
execute_on = 'initial'
timestep = 'LATEST'
[../]
[]

[AuxKernels]
#active=''
[./en]
type = SolutionAux
solution = soln
variable = P_aux
scale_factor = 1.0
execute_on = 'initial'
[../]
[]

[MultiApps]
[./sub]
type = TransientMultiApp
app_type = parrotApp
execute_on = timestep_end
input_files = plotLine_sub_P.i
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
[./tosub]
type = MultiAppProjectionTransfer
direction = to_multiapp
multi_app = sub
source_variable = P_aux
variable = u
[../]
[]
 
 
[Outputs]
 execute_on = 'timestep_end'
 [./out]
 type = Exodus
 execute_on = 'timestep_end'
 [../]
[]



