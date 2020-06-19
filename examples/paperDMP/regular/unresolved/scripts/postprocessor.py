from paraview.simple import *
import os, sys

inputName=sys.argv[1]
case=inputName[12:-2]
print(case)

pressureLineAName='./pres_line_A_'+case+'.csv'
pressureLineBName='./pres_line_B_'+case+'.csv'
pressureLineCName='./pres_line_C_'+case+'.csv'

Flux1LineAName='./Flux1_line_A_'+case+'.csv'
Flux2LineAName='./Flux2_line_A_'+case+'.csv'

reader = ExodusIIReader(FileName=inputName)
#Show(ExodusIIReader(FileName=inputName))
#reader = GetActiveSource()
#tsteps = reader.TimestepValues
#print(tsteps)

reader.GenerateObjectIdCellArray = 0
reader.GenerateGlobalElementIdArray = 0
reader.GenerateGlobalNodeIdArray = 0
reader.GlobalVariables = []

#view = GetActiveView()
#view.ViewTime = tsteps[1]
#Render()

plotOverLine1 = PlotOverLine(Input=reader, Source='High Resolution Line Source')
plotOverLine1.Source.Resolution = 20001

# A line
plotOverLine1.Source.Point1 = [0.0, 0.7, 0.0]
plotOverLine1.Source.Point2 = [1.0, 0.7, 0.0]

reader.PointVariables = ['flux_1']
SaveData(Flux1LineAName, proxy=plotOverLine1, Precision=12, WriteTimeSteps=1)

reader.PointVariables = ['flux_2']
SaveData(Flux2LineAName, proxy=plotOverLine1, Precision=12, WriteTimeSteps=1)

reader.PointVariables = ['pressure']
SaveData(pressureLineAName, proxy=plotOverLine1, Precision=12, WriteTimeSteps=0)

# B line
plotOverLine1.Source.Point1 = [0.5, 0.0, 0.0]
plotOverLine1.Source.Point2 = [0.5, 1.0, 0.0]

reader.PointVariables = ['pressure']
SaveData(pressureLineBName, proxy=plotOverLine1, Precision=12, WriteTimeSteps=0)

# C line
plotOverLine1.Source.Point1 = [0.0, 0.1, 0.0]
plotOverLine1.Source.Point2 = [0.9, 1.0, 0.0]

reader.PointVariables = ['pressure']
SaveData(pressureLineCName, proxy=plotOverLine1, Precision=12, WriteTimeSteps=0)
