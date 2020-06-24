[Problem]
type = FEProblem
solve = false
[]

[Mesh]
 file = single_${tipeMesh}_${origRef}.e
 block_id = '2 4 5 6 7'
 boundary_id = '1 2'
 boundary_name = 'inflow outflow'
 uniform_refine = ${unifSteps}
[]

[MeshModifiers]
active=''
[./subdomains_0]
 type = SubdomainBoundingBox
 bottom_left = '-51.00 -40.00001 0.0'
 top_right =   '0.0000 51.0001 0.0'
 block_id = 0
 block_name = up_block_0
 [../]
 
[./subdomains_1]
 type = SubdomainBoundingBox
 bottom_left = '-51.00 -51.0 0.0'
 top_right =   ' 0.00 -40.0 0.0'
 block_id = 2
 block_name = bottom_block_1
  execute_on=timestep_end
 [../]
 
[./subdomains_2]
 type = SubdomainBoundingBox
 bottom_left = '0.00001 -51.0 0.0'
 top_right =   '51.00000 -40.0 0.0'
 block_id = 1
 block_name = bottom_block_2
execute_on=timestep_end
 [../]
 
 
[./subdomains_3]
 type = SubdomainBoundingBox
 bottom_left = '0.00001 -40.00001 0.0'
 top_right =   '51.000 51.0001 0.0'
 block_id = 3
 block_name = up_block_3
 execute_on=timestep_end
 [../]
 
[./center_side_set_1]
 type = SideSetsBetweenSubdomains
 master_block = up_block_0
 paired_block = up_block_3
 new_boundary = 'new_side_set_1'
 depends_on = 'subdomains_3 subdomains_0'
  execute_on=timestep_end
 [../]
 
[./center_side_set_2]
 type = SideSetsBetweenSubdomains
 master_block = bottom_block_1
 paired_block = bottom_block_2
 new_boundary = 'new_side_set_2'
 depends_on = 'subdomains_1 subdomains_2'
 [../]
[]
[Variables]
[./pressure]  order = FIRST  family=LAGRANGE [../]
[]

[Kernels]
[./StressDivergenceParrot_real_x] type = Reaction variable = pressure [../]
[]

[Executioner]
 type=Steady
 solve_type=LINEAR
 line_search = 'none'
 nl_abs_tol = 1e-8
 []

[Outputs]
 file_base = refinedMesh_${unifSteps}
 exodus = true
 perf_graph = true
 xdr = true
[]

[Adaptivity]
 marker = simplemark
 steps = ${adapSteps}
 [./Markers]
 [./simplemark]
 type = BlockMarker
 blockNum = 4
 [../]
 [../]
[]
