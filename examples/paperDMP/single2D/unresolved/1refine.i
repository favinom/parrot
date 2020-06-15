[Problem]
type = FEProblem
solve = false
[]

[Mesh]
  file = refinedMesh_0000_mesh.xdr
  boundary_id = '11 22'
  boundary_name = 'inflow outflow'
  parallel_type=distributed
# uniform_refine=3
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

[./my]
type = FractureRefinement
fractureMeshModifier = fractureUserObject
refinements='6 0'
#outputFileName = mesh_6_0.e
doBoundaryRefinement = false
[../]


[]


[Executioner]
type=Steady
solve_type=LINEAR
line_search = 'none'
nl_abs_tol = 1e-8
[]

[Outputs]
 xda     = true
 perf_graph = true
 []
