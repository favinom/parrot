from paraview.simple import *
import os, sys

inputName=sys.argv[1]
case=inputName[12:-2]

CMLineAName='./CM_line_0.5'+case+'.csv'
CMLineBName='./CM_line_0.75'+case+'.csv'

reader = ExodusIIReader(FileName=inputName)

reader.GenerateObjectIdCellArray = 0
reader.GenerateGlobalElementIdArray = 0
reader.GenerateGlobalNodeIdArray = 0
reader.GlobalVariables = []

tsteps = reader.TimestepValues
print(tsteps)
plotOverLine1 = PlotOverLine(Input=reader, Source='High Resolution Line Source')

# A line
plotOverLine1.Source.Point1 = [0.0, 0.5, 0.0]
plotOverLine1.Source.Point2 = [1.0, 0.5, 0.0]

reader.PointVariables = ['CM']

SaveData(CMLineAName, proxy=plotOverLine1, Precision=12, WriteTimeSteps=1)

# B line
plotOverLine1.Source.Point1 = [0.0, 0.75, 0.0]
plotOverLine1.Source.Point2 = [1.0, 0.75, 0.0]

reader.PointVariables = ['CM']
SaveData(CMLineBName, proxy=plotOverLine1, Precision=12, WriteTimeSteps=1)

