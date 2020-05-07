[Problem]
type = FEProblem
solve = false
skip_nl_system_check = true
[]

[Mesh]
file = bcMesh.xdr
parallel_type = distributed
[]


#[Mesh]
#type = GeneratedMesh
#dim = 3
#xmin = 0.0
#xmax = 1.0
#ymin = 0.0
#ymax = 2.25
#zmin = 0.0
#zmax = 1.0
#nx = 9 # 144 # ${nx}
#ny = 20 # 324 # ${ny}
#nz = 9 # 144 # ${nx}
#parallel_type = distributed
#[]


[MeshModifiers]

[./fractureUserObject]
type = FractureUserObject
fn = 8
fx_string = '0.5,0.5,0.5,0.5,0.2,0.2,0.77,0.83'
fy_string = '1.125,0.175,1.6,1.6,2.05,2.05,2.05,2.05'
fz_string = '0.5,0.5,0.675,0.31,0.5,0.5,0.5,0.5'
fd1_string = '0.9,0.9,0.9,0.9,0.4,0.4,0.4,0.4'
fd2_string = '1.75,0.25,1.25,1.2472,0.30594,0.30594,0.3,0.3'
fd3_string = '0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01'
fa1_string = '0,0,0,0,78.6901,-78.6901,0,0'
fa2_string = '0,90,0,0,-90,-90,-90,-90'
fa3_string = '0,0,16.2602,-15.8192,90,-90,0,0'
[../]

[]


[Executioner]
 type=Steady
 solve_type=LINEAR
 line_search = 'none'
 nl_abs_tol = 1e-8
 []

[UserObjects]
[./soln]
type = RefineMesh
fractureMeshModifier = fractureUserObject
refinements = '1 1 1 1 2 0'
filename = 'adaptedMesh.e'
flag=2
[../]
[]