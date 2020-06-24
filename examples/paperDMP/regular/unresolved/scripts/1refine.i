[Problem]
type = FEProblem
solve = false
[]

[Mesh]
file = refinedMesh_0000_mesh.xdr
parallel_type=distributed
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

[./my]
type = FractureRefinement
fractureMeshModifier = fractureUserObject
refinements='${mRefLev} 0'
# outputFileName = mesh_${mResName}_${mRefLevName}_${mUmr}.e
doBoundaryRefinement = false
[../]


[]

[Outputs]
 file_base = scripts/1refinedMesh
 #exodus = true
 perf_graph = true
 xdr = true
[]
 
[Executioner]
type=Steady
solve_type=LINEAR
line_search = 'none'
nl_abs_tol = 1e-8
[]

