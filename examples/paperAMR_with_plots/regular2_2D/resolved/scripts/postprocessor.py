from paraview.simple import *
import os, sys

inputName=sys.argv[1]
case=inputName[12:-2]

pressureLineAName='./pres_line_A'+case+'.csv'
pressureLineBName='./pres_line_B'+case+'.csv'
pressureLineCName='./pres_line_C'+case+'.csv'

CMLineAName='./CM_line_A'+case+'.csv'
CMLineBName='./CM_line_B'+case+'.csv'
CMLineCName='./CM_line_C'+case+'.csv'

reader = ExodusIIReader(FileName=inputName)

reader.GenerateObjectIdCellArray = 0
reader.GenerateGlobalElementIdArray = 0
reader.GenerateGlobalNodeIdArray = 0
reader.GlobalVariables = []

plotOverLine1 = PlotOverLine(Input=reader, Source='High Resolution Line Source')

# A line
plotOverLine1.Source.Point1 = [0.5, 0.0, 0.0]
plotOverLine1.Source.Point2 = [1.5, 1.0, 0.0]

reader.PointVariables = ['CM']

SaveData(CMLineAName, proxy=plotOverLine1, Precision=12, WriteTimeSteps=1)

reader.PointVariables = ['pressure']
SaveData(pressureLineAName, proxy=plotOverLine1, Precision=12, WriteTimeSteps=0)

# B line
plotOverLine1.Source.Point1 = [0.0, 0.75, 0.0]
plotOverLine1.Source.Point2 = [1.0, 0.75, 0.0]

reader.PointVariables = ['CM']
SaveData(CMLineBName, proxy=plotOverLine1, Precision=12, WriteTimeSteps=1)


reader.PointVariables = ['pressure']
SaveData(pressureLineBName, proxy=plotOverLine1, Precision=12, WriteTimeSteps=0)

# C line
plotOverLine1.Source.Point1 = [0.0, 0.5, 0.0]
plotOverLine1.Source.Point2 = [1.0, 0.5, 0.0]

reader.PointVariables = ['CM']
SaveData(CMLineCName, proxy=plotOverLine1, Precision=12, WriteTimeSteps=1)


reader.PointVariables = ['pressure']
SaveData(pressureLineCName, proxy=plotOverLine1, Precision=12, WriteTimeSteps=0)
