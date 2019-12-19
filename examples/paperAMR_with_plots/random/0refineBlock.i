[Problem]
type = FEProblem
solve = false
[]

[Mesh]
 type = GeneratedMesh
 xmin=-200.0
 xmax= 200.0
 ymin=-200.0
 ymax= 200.0
 nx = 32
 ny = 32
 dim = 2
[]

[MeshModifiers]

[./fractureUserObject]
type = FractureUserObject
fn = 62
fx_string = '193.18,37.424,-127.36,-156.83,109.06,-187.32,-85.269,184.06,40.5,-13.475,-71.621,-95.027,156.28,160.95,81.93,13.113,157.8,170.28,123.92,-20.976,-183.8,-198.56,-111.82,-30.973,-36.042,115.02,-72.734,88.463,120.44,128.72,50.982,-163.23,25.064,-135.44,-26.407,83.367,167.71,-117.29,-106.23,-65.607,18.594,78.298,-169.8,-120.06,-152.33,180,-121,106.92,131.5,-159.38,1.5643,75.768,25.927,-114.52,-82.963,-45.737,-104.98,164.82,175.69,115.32,146.23,39.169
83.855,153.45,76.575,-128.44,97.627,56.936,82.477,177.25,-71.525,180.61,12.659,5.6952,-3.7316,-139.11,-47.603,-58.14,-114.33,159.63,98.093,55.184,-30.191,-86.845,143.16,177.49,-51.525,187.77,-16.034,-170.74,-80.147,26.272,85.565,197.21,100.35,54.182,113.93,-182.72,-151.5,42.946,124.24,-142.06,20.331,119.93,155.21,-134.81,182.11,49.873,144.8,25.066,-2.4792,67.985,12.37,83.149,-180.34,-185.45,132.97,142.78,-26.286,-165.43,-183.9,-21.086,-25.85,164.61
-129.83,-30.21,-114.49,98.621,-22.838,26.748,14.082,-154.92,-86.283,186,149.28,-21.29,122.33,-122.97,-169.44,109.25,-198.46,17.933,-65.143,83.778,16.086,55.44,161.89,-199.47,-109.26,-127.81,-105.64,37.945,110.25,-178.25,-77.442,-161.43,-196.31,137.58,10.094,-141.62,36.778,-60.948,-19.49,-102.05,77.815,64.242,24.979,-185.25,-23.64,-117.22,-149.74,158.86,36.311,-61.041,-48.351,62.769,148.91,-22.182,28.21,69.93,-189.98,7.9804,-134.54,-160.15,-131.17,-146.81
143.32,42.903,-80.758,-180.21,-178.8,-49.436,-122.72,59.315,-25.874,106.11,-178.18,-63.516,-69.424,116.39,-35.664,152.67,152.23,160.45,33.73,197.05,181.53,36.83,-83.192,192.51,-21.588,172.24,-188.9,144.79,21.108,-95.994,-94.525,-74.741,-9.2861,112.91,-67.489,-106.67,-56.139,87.095,-100.01,-48.411,171,7.9746,-122.26,-90.91,150.3,-156.14,58.232,-52.768,-182.04,67.798,107.56,111.59,144.01,-69.368,-171.8,43.239,-92.541,21.052,-24.751,122.33,-198.28,-183.13
163.76,-171.69,107.33,-171.49,-164.87,-114.98,75.774,-7.6784,161.5,29.813,0.16031,135.66,19.952,-175.72,-142.8,93.637,-105.95,-179.27,-12.419,172.88,-116.44,-69.884,90.348,28.078,-93.511,-181.94,63.381,-20.48,21.877,35.641,166.4,114.15,-99.868,-94.173,-27.363,-101.32,87.723,-188.8,182.18,-91.872,177.68,-67.305,-111.42,-107.95,145.98,26.34,-25.073,-64.908,25.783,-81.879,-86.995,32.188,-64.697,-85.082,-93.475,-24.23,192.06,74.06,-69.733,63.782,60.999,-39.584'
fy_string = '-67.0608,112.8069,78.2961,27.8333,28.6952,136.3192,147.5871,117.9536,91.1739,48.1563,-6.0869,110.5975,-0.33828,-112.2095,188.1677,-178.9852,72.3103,133.4351,-114.7376,-8.5073,46.8002,21.6364,46.2736,-175.1089,165.1034,-37.3441,51.2528,150.3301,139.2483,-6.2823,10.5859,51.719,-3.5615,61.1705,-117.5612,104.1047,84.9939,38.9261,-72.3578,31.2006,63.7639,-58.7467,112.0272,180.016,-9.2725,-15.4259,9.7291,-18.9361,-119.2662,-123.561,184.7674,-160.7251,-55.6731,198.7802,-157.9315,187.4522,-173.9624,187.4522,-173.9624,-197.9318,151.9745,-89.6609'
fa1_string = '131.1072,147.675,41.1714,88.4972,88.4972,96.6395,141.7574,79.0334,39.0899,133.5214,133.5214,103.3064,128.1842,107.6244,135.3764,135.3764,132.8385,81.0055,113.5286,73.9205,73.1704,119.2585,119.2585,136.4225,136.4225,57.2585,61.3704,99.1981,99.1981,60.8217,60.8217,149.695,102.4057,125.3422,73.0231,107.7189,107.7189,48.2934,117.9998,117.3385,53.0861,103.7523,34.4829,96.8455,139.5299,50.1438,50.1438,113.8547,48.1872,71.4659,34.7866,34.7866,102.5925,124.6825,124.6825,36.0741,36.0741,36.0741,36.0741,47.0323,47.0323,131.491'
fd1_string = '150.2898,139.245,141.012,18.2782,47.4513,173.0121,113.9766,135.9827,25.1432,44.9343,112.606,59.1174,30.3359,46.8148,33.249,59.0522,74.4429,178.9839,164.7223,175.733,24.4857,26.4948,74.3223,68.718,96.3406,136.7864,51.2114,1.8481,136.806,24.5226,44.6763,158.7088,154.3843,52.245,11.3319,46.9829,78.6029,32.177,163.0294,145.8514,148.9826,50.3868,148.3879,27.3985,50.5654,7.5405,70.9628,151.2911,76.761,155.7068,37.0947,95.6427,141.3407,4.2874,147.8606,31.049,64.4291,31.049,64.4291,6.0687,140.9227,131.0122'
fd2_string = '0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5'
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
 file_base = refinedMesh
# exodus = true
 perf_graph = true
 xdr = true
[]

[Adaptivity]
 marker = simplemark
 steps = ${adapSteps}
 [./Markers]
 [./simplemark]
 type = RegionMarker
 regionMeshModifier = fractureUserObject
 [../]
 [../]
[]
