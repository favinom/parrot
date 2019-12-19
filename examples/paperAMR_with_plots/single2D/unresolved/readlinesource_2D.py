# state file generated using paraview version 5.4.1

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------
import os, sys

# folder = sys.argv[-1]
folder = os.getcwd()
masterfiles = []
slavefiles = []
name = input("Enter file name")
print(name)
for file in os.listdir(folder):
    if file.find(name)>-1:
        masterfiles.append(os.path.join(folder, file))

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

##################
# read master info
##################

# create a new 'ExodusIIReader'
mastere1 = ExodusIIReader(FileName=sorted(masterfiles))

writer = CreateWriter("%s/results_line_1_cm.csv"%folder, mastere1)
writer.UpdatePipeline()


