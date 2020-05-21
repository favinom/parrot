[GlobalParams]
conservative = false
[]

[Problem]
type = ParrotProblem3
use_AFC = true
operator_userobject = storeOperatorsUO
solver_type = 1
#antidiffusive_fluxes=antidiffusive_fluxes
[]

#[Problem]
#type = ParrotProblem
#use_AFC = true
#operator_userobject = storeOperatorsUO
#antidiffusive_fluxes = antidiffusive_fluxes
#[]


[Mesh]
 type = GeneratedMesh
 xmin= 0.0
 xmax= 700.0
 ymin= 0.0
 ymax= 600.0
 nx = 70
 ny = 60
 dim = 2
 parallel_type = distributed
[]

[MeshModifiers]
[./fractureUserObject]
type = FractureUserObject
fn = 63
fx_string =  '313.2676,260.8653,265.1721,269.8985,346.7056,396.9458,186.6971,152.8704,62.9313,413.3329,460.6934,382.6412,364.9141,385.0224,354.0661,534.7769,540.8623,580.3708,422.2058,377.4094,441.2443,537.5565,613.7643,220.6222,180.3054,166.8665,113.1108,204.9435,445.724,496.1199,493.88,527.4774,50.3959,415.1575,568.5499,558.8264,565.3749,493.7388,500.6841,524.4967,500.8826,601.0937,615.3813,629.2719,403.648,407.6168,392.5355,481.6341,477.2685,474.8872,365.7464,342.7276,326.4557,384.2011,356.0229,448.2966,517.1544,513.1858,563.1921,570.7327,577.2812,578.8687,586.2108'
fy_string =  '231.0968,492.9808,503.2714,336.242,474.2399,469.0665,410.1637,430.8623,407.776,154.277,122.527,138.9312,206.4,227.5667,303.5023,366.2087,357.2129,337.3153,487.1608,492.7603,425.5657,442.3644,536.5078,166.8666,469.2422,370.6902,453.5635,435.6449,212.7829,193.7444,313.5748,207.1834,244.1404,397.9572,378.7087,353.7055,337.4337,387.8369,428.5166,370.3743,361.643,513.0512,504.5184,489.4371,573.5747,566.4309,576.9481,493.2074,207.0202,190.3515,161.1811,147.6873,165.5467,255.0422,232.2218,183.0092,211.7828,191.3437,186.3827,188.5656,194.717,305.247,303.2626'
fa1_string = '61.088,-62.7136,-61.0239,89.7678,-85.8502,77.1063,-86.4767,89.4608,72.6581,77.2987,89.2631,-49.109,-51.3966,-52.4521,-51.8817,51.9682,53.8418,42.5107,46.8677,19.0256,-31.5042,-32.981,-36.3629,-37.1019,37.1763,-28.8278,-31.3995,-30.1317,-55.0975,-59.6764,-22.4022,-59.0362,86.9335,-38.6598,39.5749,52.0069,52.5641,56.7863,71.565,53.3625,53.8108,28.6106,26.5651,25.8769,45.9547,51.2258,55.1094,41.4235,40.3177,42.3575,52.125,53.6731,58.3176,57.5289,-48.7376,-53.0412,-56.7437,-56.352,-82.7468,-83.5169,-83.4802,42.8388,44.1696'
fd1_string = '180.5979,49.5315,28.1273,376.9243,252.1814,268.6406,375.5885,338.2904,224.2987,247.8957,164.5848,261.0992,167.9288,178.8712,43.7195,135.7064,119.2837,116.5242,291.5725,68.708,162.8766,250.9897,214.1761,553.24,418.8521,329.809,202.0559,370.3363,234.8716,137.5256,317.3675,143.663,125.6099,83.861,50.4595,101.2229,104.4628,145.6324,70.2818,137.0027,113.5939,59.6741,64.7833,59.1079,33.6807,31.0528,36.7726,35.9911,147.8261,91.3031,171.9529,66.9956,119.393,67.2701,55.9651,204.6294,159.9415,149.7,66.0128,52.7247,41.9433,59.5366,58.0971'
fd2_string = '0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01'
[../]


[./my]
type = FractureRefinement
fractureMeshModifier = fractureUserObject
refinements='${mRefLev} ${mUmr}'
# outputFileName= 'ciao.e'
doBoundaryRefinement = true
[../]
[]

[Variables]
[./CM] [../]
[]
 
[AuxVariables]
[./pressure] [../]
[./correction] [../]
[]

[Materials]
[./conductivity1] type = FlowAndTransport fractureMeshModifier =  fractureUserObject
phi = ${phiIn} phiFrac = ${phiIn}
k = 8.64e-10 kFrac = 8.64e-4
pressure = pressure
[../]
[]

[Kernels]
active='time upwind'
[upwind] type = Advection variable = CM [../]
[./time] type = PorosityTimeDerivative variable = CM lumping = true [../]
[]

[BCs]
[./u_injection_left] type = DirichletBC boundary = left variable = CM value='1.0' [../]
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

 dt = 73
 end_time = 36500

petsc_options_iname=' -ksp_type            '   # -mat_view
petsc_options_value='  ksp_parrot_preonly  '   # ::ascii_matlab


[./Quadrature] order= NINTH type = GRID [../]

[]

[Outputs]
 file_base = AdvectionOut_${mRefLevName}_${mUmr}
[./exod]
 type  = Exodus
 interval = 5
[../]
[./exod2]
 type = CSV
[../]

[]


[UserObjects]
[./soln]
type = SolveDiffusion2
execute_on = 'initial'
block_id='0'
value_p ='0.001 1000'
boundary_D_bc='3 1'
value_D_bc='1013250.0 0.0'
boundary_N_bc=' '
value_N_bc=' '
aux_variable=pressure
fractureMeshModifier = fractureUserObject
output_file=DiffusionOut2_${mRefLevName}_${mUmr}.e
solver_type = 1
[../]
[./storeOperatorsUO]
type = StoreOperators
[../]
[./MassAssembly]
type = AssembleMassMatrix
operator_userobject = storeOperatorsUO 
block_id = '0'
value_p = '${phiIn} ${phiIn}'
execute_on = 'initial'
constrain_matrix = true
fractureMeshModifier = fractureUserObject

dc_boundaries = '3'
value_D_bc='1'
dc_variables='CM'
[../]

[./assembleVolumeVectors]
type=AssembleVolumeVectors
RegionMeshModifier = fractureUserObject
execute_on = 'initial'
[../]


#[./antidiffusive_fluxes]
# type = AntidiffusiveFluxes
# operator_userobject = storeOperatorsUO
# execute_on = 'timestep_end'
# dc_boundaries = '3 1'
# WriteCorrection=true
#[../]

[]



[Postprocessors]
[./int0] type = IntegralSolutionOverRegionFast region = 0 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int1] type = IntegralSolutionOverRegionFast region = 1 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int2] type = IntegralSolutionOverRegionFast region = 2 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int3] type = IntegralSolutionOverRegionFast region = 3 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int4] type = IntegralSolutionOverRegionFast region = 4 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int5] type = IntegralSolutionOverRegionFast region = 5 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int6] type = IntegralSolutionOverRegionFast region = 6 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int7] type = IntegralSolutionOverRegionFast region = 7 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int8] type = IntegralSolutionOverRegionFast region = 8 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int9] type = IntegralSolutionOverRegionFast region = 9 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int10] type = IntegralSolutionOverRegionFast region = 10 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int11] type = IntegralSolutionOverRegionFast region = 11 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int12] type = IntegralSolutionOverRegionFast region = 12 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int13] type = IntegralSolutionOverRegionFast region = 13 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int14] type = IntegralSolutionOverRegionFast region = 14 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int15] type = IntegralSolutionOverRegionFast region = 15 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int16] type = IntegralSolutionOverRegionFast region = 16 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int17] type = IntegralSolutionOverRegionFast region = 17 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int18] type = IntegralSolutionOverRegionFast region = 18 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int19] type = IntegralSolutionOverRegionFast region = 19 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int20] type = IntegralSolutionOverRegionFast region = 20 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int21] type = IntegralSolutionOverRegionFast region = 21 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int22] type = IntegralSolutionOverRegionFast region = 22 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int23] type = IntegralSolutionOverRegionFast region = 23 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int24] type = IntegralSolutionOverRegionFast region = 24 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int25] type = IntegralSolutionOverRegionFast region = 25 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int26] type = IntegralSolutionOverRegionFast region = 26 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int27] type = IntegralSolutionOverRegionFast region = 27 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int28] type = IntegralSolutionOverRegionFast region = 28 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int29] type = IntegralSolutionOverRegionFast region = 29 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int30] type = IntegralSolutionOverRegionFast region = 30 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int31] type = IntegralSolutionOverRegionFast region = 31 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int32] type = IntegralSolutionOverRegionFast region = 32 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int33] type = IntegralSolutionOverRegionFast region = 33 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int34] type = IntegralSolutionOverRegionFast region = 34 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int35] type = IntegralSolutionOverRegionFast region = 35 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int36] type = IntegralSolutionOverRegionFast region = 36 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int37] type = IntegralSolutionOverRegionFast region = 37 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int38] type = IntegralSolutionOverRegionFast region = 38 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int39] type = IntegralSolutionOverRegionFast region = 39 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int40] type = IntegralSolutionOverRegionFast region = 40 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int41] type = IntegralSolutionOverRegionFast region = 41 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int42] type = IntegralSolutionOverRegionFast region = 42 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int43] type = IntegralSolutionOverRegionFast region = 43 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int44] type = IntegralSolutionOverRegionFast region = 44 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int45] type = IntegralSolutionOverRegionFast region = 45 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int46] type = IntegralSolutionOverRegionFast region = 46 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int47] type = IntegralSolutionOverRegionFast region = 47 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int48] type = IntegralSolutionOverRegionFast region = 48 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int49] type = IntegralSolutionOverRegionFast region = 49 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int50] type = IntegralSolutionOverRegionFast region = 50 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int51] type = IntegralSolutionOverRegionFast region = 51 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int52] type = IntegralSolutionOverRegionFast region = 52 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int53] type = IntegralSolutionOverRegionFast region = 53 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int54] type = IntegralSolutionOverRegionFast region = 54 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int55] type = IntegralSolutionOverRegionFast region = 55 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int56] type = IntegralSolutionOverRegionFast region = 56 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int57] type = IntegralSolutionOverRegionFast region = 57 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int58] type = IntegralSolutionOverRegionFast region = 58 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int59] type = IntegralSolutionOverRegionFast region = 59 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int60] type = IntegralSolutionOverRegionFast region = 60 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int61] type = IntegralSolutionOverRegionFast region = 61 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]
[./int62] type = IntegralSolutionOverRegionFast region = 62 doDomainSize = 0 VolumeUserObject = assembleVolumeVectors [../]

[./reg0] type = IntegralSolutionOverRegionFast region = 0 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg1] type = IntegralSolutionOverRegionFast region = 1 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg2] type = IntegralSolutionOverRegionFast region = 2 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg3] type = IntegralSolutionOverRegionFast region = 3 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg4] type = IntegralSolutionOverRegionFast region = 4 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg5] type = IntegralSolutionOverRegionFast region = 5 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg6] type = IntegralSolutionOverRegionFast region = 6 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg7] type = IntegralSolutionOverRegionFast region = 7 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg8] type = IntegralSolutionOverRegionFast region = 8 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg9] type = IntegralSolutionOverRegionFast region = 9 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg10] type = IntegralSolutionOverRegionFast region = 10 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg11] type = IntegralSolutionOverRegionFast region = 11 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg12] type = IntegralSolutionOverRegionFast region = 12 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg13] type = IntegralSolutionOverRegionFast region = 13 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg14] type = IntegralSolutionOverRegionFast region = 14 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg15] type = IntegralSolutionOverRegionFast region = 15 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg16] type = IntegralSolutionOverRegionFast region = 16 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg17] type = IntegralSolutionOverRegionFast region = 17 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg18] type = IntegralSolutionOverRegionFast region = 18 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg19] type = IntegralSolutionOverRegionFast region = 19 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg20] type = IntegralSolutionOverRegionFast region = 20 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg21] type = IntegralSolutionOverRegionFast region = 21 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg22] type = IntegralSolutionOverRegionFast region = 22 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg23] type = IntegralSolutionOverRegionFast region = 23 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg24] type = IntegralSolutionOverRegionFast region = 24 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg25] type = IntegralSolutionOverRegionFast region = 25 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg26] type = IntegralSolutionOverRegionFast region = 26 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg27] type = IntegralSolutionOverRegionFast region = 27 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg28] type = IntegralSolutionOverRegionFast region = 28 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg29] type = IntegralSolutionOverRegionFast region = 29 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg30] type = IntegralSolutionOverRegionFast region = 30 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg31] type = IntegralSolutionOverRegionFast region = 31 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg32] type = IntegralSolutionOverRegionFast region = 32 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg33] type = IntegralSolutionOverRegionFast region = 33 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg34] type = IntegralSolutionOverRegionFast region = 34 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg35] type = IntegralSolutionOverRegionFast region = 35 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg36] type = IntegralSolutionOverRegionFast region = 36 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg37] type = IntegralSolutionOverRegionFast region = 37 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg38] type = IntegralSolutionOverRegionFast region = 38 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg39] type = IntegralSolutionOverRegionFast region = 39 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg40] type = IntegralSolutionOverRegionFast region = 40 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg41] type = IntegralSolutionOverRegionFast region = 41 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg42] type = IntegralSolutionOverRegionFast region = 42 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg43] type = IntegralSolutionOverRegionFast region = 43 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg44] type = IntegralSolutionOverRegionFast region = 44 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg45] type = IntegralSolutionOverRegionFast region = 45 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg46] type = IntegralSolutionOverRegionFast region = 46 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg47] type = IntegralSolutionOverRegionFast region = 47 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg48] type = IntegralSolutionOverRegionFast region = 48 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg49] type = IntegralSolutionOverRegionFast region = 49 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg50] type = IntegralSolutionOverRegionFast region = 50 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg51] type = IntegralSolutionOverRegionFast region = 51 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg52] type = IntegralSolutionOverRegionFast region = 52 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg53] type = IntegralSolutionOverRegionFast region = 53 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg54] type = IntegralSolutionOverRegionFast region = 54 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg55] type = IntegralSolutionOverRegionFast region = 55 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg56] type = IntegralSolutionOverRegionFast region = 56 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg57] type = IntegralSolutionOverRegionFast region = 57 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg58] type = IntegralSolutionOverRegionFast region = 58 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg59] type = IntegralSolutionOverRegionFast region = 59 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg60] type = IntegralSolutionOverRegionFast region = 60 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg61] type = IntegralSolutionOverRegionFast region = 61 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[./reg62] type = IntegralSolutionOverRegionFast region = 62 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]

[]