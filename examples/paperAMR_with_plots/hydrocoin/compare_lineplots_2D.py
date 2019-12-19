import os, sys
import json
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from matplotlib.ticker import FormatStrFormatter
import numpy as np

def plotComparison(files, colord, outname, zoom=False):
    fig = plt.figure(num=None, figsize=(5, 4), dpi=1200, facecolor='w', edgecolor='k')
    plt.tick_params(axis='both', which='major')
    
    for file in files:
        filename = file[0]
        lines = open(filename,'r').readlines()
        dataDict = {header:i for i, header in enumerate(lines[0][:-1].split(','))}
        dataDict['values'] = np.array([[float(item) for item in line.split(',')] for line in lines[1:]])
        #dataDict2 = {header:i for i, header in enumerate(lines[0][:-1].split(','))}
        #dataDict2['values'] = np.array([[float(item) for item in line.split(',')] for line in lines[5:]])
    
    
        Y =  dataDict['values'][:,dataDict[file[1]]]
        X =  dataDict['values'][:,dataDict[file[2]]]

    
        labelname = file[3]
        if labelname.find('Reference') >-1:
            marks = int(round(len(X)/30.))
            if zoom:
                lines = plt.plot(X, Y, label=labelname, marker='x', markersize = 20, markevery=marks, markerfacecolor='None', color='k')
            else:
                lines = plt.plot(X, Y, label=labelname, marker='*', markersize = 7, markevery=marks, markerfacecolor='k', color='k')
# lines1 = plt.plot(X1, Y1, label=labelname, marker='o', markersize = 20, markevery=marks, markerfacecolor='None', color='k')
        elif filename.find('3D') >-1:
            lines = plt.plot(X, Y, label=labelname, color='r', alpha=colord[file[4]])
        else:
            lines = plt.plot(X, Y, label=labelname, color='b', alpha=colord[file[4]])
        if len(file)>5:
            lines[0].set_color(file[5])

    if not zoom:
        plt.xlabel('Arc length [m]')
        plt.ylabel('Pressure [m]')
#        plt.ylim([106, 120])

    if zoom:
        outname = outname + '_zoom'
        plt.xlim([.6, 1.])
        plt.ylim([1.132, 1.148])
        plt.locator_params(axis='y', nbins=6)
    else:
        #plt.legend(loc = 9, bbox_to_anchor=(.50, 1.28), frameon=True, ncol=3)
        plt.legend(loc = 8, bbox_to_anchor=(.50, -0.42), frameon=True, ncol=3)
        #plt.legend(frameon=True, ncol=2, loc='upper right')
        if outname.find('y')>-1:
            line_name = 'AB'
        if outname.find('x')>-1:
            line_name = 'BB$^{\prime}$'
        plt.title('Line %s'%line_name)
    plt.grid(True)
    
    plt.subplots_adjust(left=0.20, right=0.95, top=0.90, bottom=0.15)
    plt.savefig('%s_2D.eps'%outname, bbox_inches='tight')


reference_base = '/Users/mariagiuseppinanestola/Documents/ICS/Benchmarking/results_hydrocoin/ottimo/fracture-flow-master/hydrocoin/results'
reference_d = os.path.join(reference_base, 'mfd')
ref_y = 'mfd_hydrocoin_200.csv'

working_d = os.getcwd()


# solvers = [
# #           'monolithic_quad4',
#            'monolithic',
# #           'p0',
# #           'static'
#            ]

out_y = 'results_y_eq_-200_2D.csv'

colord = {257: 1.0,
          129: 0.7,
          65:  0.5,
          33:  0.3,
         }
#
#files = [
#        (os.path.join(reference_d, ref_y), '"pressure"', '"arc_length"', 'Reference'),
#        (os.path.join(working_d, meshes[3], solvers[0], out_y), '"u_m"', '"arc_length"', '$\sqrt{\# cells}$=257', 257),
#        (os.path.join(working_d, meshes[2], solvers[0], out_y), '"u_m"', '"arc_length"', '$\sqrt{\# cells}$=129', 129),
#        (os.path.join(working_d, meshes[1], solvers[0], out_y), '"u_m"', '"arc_length"', '$\sqrt{\# cells}$=65', 65),
#        (os.path.join(working_d, meshes[0], solvers[0], out_y), '"u_m"', '"arc_length"', '$\sqrt{\# cells}$=33', 33),
#        ]
#
#plotComparison(files, colord, 'comparison_y_0.7')


#f = open("colormap.json")
#data = json.load(f)
#f.close()
#colors = data[0][u'RGBPoints']
#diff = 4
#colm = []
#for i in range(len(colors)//diff):
#    colm.append(tuple(colors[i+1:i+diff]))

cm = plt.get_cmap('Set1')
NUM_COLORS = 9
colors = lambda j : cm(j//1*float(1)/NUM_COLORS)
colorsidx = [0, 1, 2, 3, 4, 8, 6]
colormap = lambda j : [colors(idx) for idx in colorsidx][j]

files = [
        (os.path.join(reference_d, ref_y), '"pressure"', '"arc_length"', 'Reference'),
        ]
fileend = '_hydrocoin_200.csv'
# ### Box
# path = os.path.join(reference_base, 'boxdfm', 'boxdfm'+fileend)
# files.append((path, '"pw"', '"Points:0"', 'Box', 257, colormap(0)))
# ### TPFA
# path = os.path.join(reference_base, 'ccdfm', 'ccdfm'+fileend)
# files.append((path, '"P"', '"Points:0"', 'TPFA', 257, colormap(1)))
# ### MPFA
# path = os.path.join(reference_base, 'ccdfm', 'mpfa'+fileend)
# files.append((path, '"P"', '"Points:0"', 'MPFA', 257, colormap(2)))
# ### EDFM
# path = os.path.join(reference_base, 'edfm', 'edfm'+fileend)
# files.append((path, '"Pressure"', '"arc_length"', 'EDFM', 257, colormap(3)))
# ### Flux-Mortar
# path = os.path.join(reference_base, 'mortardfm', 'mortardfm'+fileend)
# files.append((path, '"u"', '"Points:0"', 'Flux-Mortar', 257, colormap(4)))
### P-XFEM
path = os.path.join(reference_base, 'pxfem', 'pxfem'+fileend)
#files.append((path, '"pressure"', '"arc_length"', 'P-XFEM', 257, colormap(5)))
### D-XFEM
path = os.path.join(reference_base, 'dxfem', 'dxfem'+fileend)
#files.append((path, '"Pressure"', '"arc_length"', 'D-XFEM', 257, colormap(6)))

# files.append((os.path.join(working_d, solvers[0], out_y), '"u_m"', '"arc_length"', 'LM-L$^2$', 257, 'k'))
files.append((os.path.join(working_d, out_y), '"u"', '"arc_length"', 'AMR', 257, 'r'))

plotComparison(files, colord, 'comparison_y_200')



