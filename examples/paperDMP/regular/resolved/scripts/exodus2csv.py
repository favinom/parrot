from paraview.simple import *
import os, sys

inputName=sys.argv[1]
case=inputName[12:-2]
print(case)

Flux1LineAName='./Flux2D_'+case+'_.csv'

reader = ExodusIIReader(FileName=inputName)

reader.GenerateObjectIdCellArray = 0
reader.GenerateGlobalElementIdArray = 0
reader.GenerateGlobalNodeIdArray = 0
reader.GlobalVariables = []


reader.PointVariables = ['flux_1']
SaveData(Flux1LineAName, proxy=reader, Precision=12, WriteTimeSteps=1)

