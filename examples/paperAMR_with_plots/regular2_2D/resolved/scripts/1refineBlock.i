[Problem]
type = FEProblem
solve = false
[]

[Mesh]
 file = ../adapted_${mBe}_${mFe}.xda
 parallel_type = distributed
[]

 [MeshModifiers]
 
 [./aa1]
 type = FractureUserObject
 fn = 6
 fx_string = '0.5,0.5,0.75,0.75,0.625,0.625'
 fy_string = '0.5,0.5,0.75,0.75,0.625,0.625'
 fa1_string = '0.0,90.0,0.0,90.0,0.0,90.0'
 fd1_string = '1.0,1.0,0.5,0.5,0.25,0.25'
 fd2_string = '1.0e-4,1.0e-4,1.0e-4,1.0e-4,1.0e-4,1.0e-4'
 [../]
 
 [./assign_id]
 type = RegionSubdomain
 block_id = 1
 regionMeshModifier = aa1
 [../]

[./output]
type = FractureRefinement
fractureMeshModifier = aa1
refinements='0 0' 
outputFileName = refined_${mBe}_${mFe}.xdr
# doBoundaryRefinement = true
[../]


 []

[Executioner]
type=Steady
solve_type=LINEAR
[]
