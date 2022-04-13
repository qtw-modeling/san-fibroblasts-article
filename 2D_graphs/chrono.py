import itertools as it
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
#from matplotlib_scalebar.scalebar import ScaleBar
#from scalebar import *
from matplotlib import colors as crs
from scipy.interpolate import interp1d
from scipy.ndimage.filters import gaussian_filter
from transformCmap import *
from matplotlib.font_manager import FontProperties
#from common import *

import matplotlib.dates as mdates
import matplotlib
import sys
import os
import numba as nb
from matplotlib.ticker import FormatStrFormatter

#from matplotlib import rcParams
#rcParams.update({'figure.autolayout': True})

font = {'family': 'normal', 'weight': 'normal', 'size': 30}
matplotlib.rc('font', **font)
import math as m


dt= 1e-1# 1e-4
SZ = 40
PROBstring = '00'

'''
Cmap = []
for i in range(256):
    Cmap.append(cm.autumn(i))

Cmap = np.array(Cmap)



Cmap = TransformColormap(rpd, Cmap)
Cmap = crs.LinearSegmentedColormap.from_list('my_colormap', Cmap, 256)
print Cmap
'''



dirName = '/media/karpaev/Elements/Chrono_copy2/PROB' + PROBstring + '/'

files = ['f11.dat', 'f12.dat', 'f13.dat', 'f21.dat', 'f22.dat', 'f23.dat']
#files = ['v11.dat', 'v12.dat', 'v13.dat', 'v21.dat', 'v22.dat', 'v23.dat']
titles = ['A', 'B', 'C', 'D', 'E', 'F']


fileList = [dirName + item for item in files]

#fileList.append('/home/karpaev/Videos/big_gaps/PROB00/chrono_big_no_ACh.dat')

print fileList
mapList = []

for file in fileList:
    mapList.append(np.reshape(np.fromfile(file, dtype=np.float32, count=-1, sep=""), (200,200)))

# processing nans
for k in range(200):
    for m in range(200):
        if np.isnan(mapList[3][k, m]):
            mapList[3][k, m] = 0.



#print mapList[3]

Levels=[np.arange(1492,1527,5), np.arange(1527,1577,10), np.arange(1597,1697,20), np.arange(1737,3000,40)]
Levels=list(it.chain.from_iterable(Levels))
Levels= (np.array(Levels)- 1492)*dt
#print 'Levels = ', Levels


cellLength=0.07
yMin = 0
yMax = 200
area_length = 14.
scaling_k = area_length / yMax
tickLabelsX = [0, 7, 14]
tickLabelsY = [7, 14]
tikz = np.zeros(yMax)

for i in range(200):
    if i < 200:
        tikz[i] = scaling_k * i


'''fig, axes = plt.subplots(2, 3)
fig.set_size_inches(16,16)
#fig.subplots_adjust(hspace=0.1, wspace=0.001)
axes = axes.ravel()'''

fig = plt.figure(figsize=(13.6,9))
axes = [fig.add_subplot(2,3,i+1) for i in range(6)]


Cmap = []
for i in range(len(mapList)):
    if i == 3:
        #print mapList[3]
        MIN = np.max(np.max(mapList[3])) # initial value
        for k in range(200):
            for m in range(200):
                elem = mapList[3][k,m]
                if (elem < MIN) and (elem !=0):
                    MIN = elem

        mapList[3][:,:] = (mapList[3][:,:] - MIN)*dt
        print mapList[3]

        for k in range(200):
            for m in range(200):
                if mapList[3][k,m] < 0:
                    mapList[3][k,m] = 0
    else:
        mapList[i][:,:] = (mapList[i][:,:] - np.amin(np.amin(mapList[i])))*dt

    #mapList[i] = gaussian_filter(mapList[i], sigma=0) # gaussian filter looses critical details


    CS = axes[i].contourf(mapList[i], levels=Levels, inline=1, extent=[0, 14, 0,14])
    if i == 0:
        x= CS.levels
        clim = [np.min(x), np.max(x)]
        dx = np.min(np.diff(x))
        y = np.arange(clim[0], clim[-1], dx)
        Nx = len(x)

        Cmap1 =[]
        for ii in range(Nx):
            Cmap1.append(cm.pink(ii / float(Nx - 1)))#YlGnBu_r
            #print Cmap1[i]
        #cmap = np.array(cmap)
        Cmap1 = np.array(Cmap1)
        #print Cmap1[:,0]
        #print cmap
        #print interp1d(x[:],Cmap1[:,0])(y[:])
        cmap2 = np.array([interp1d(x[:],Cmap1[:,0])(y[:]), interp1d(x[:],Cmap1[:,1])(y[:]), interp1d(x[:],Cmap1[:,2])(y[:]),\
                                                                                            interp1d(x[:],Cmap1[:,3])(y[:])])
        Cmap = crs.LinearSegmentedColormap.from_list('my_colormap', np.transpose(cmap2), 256)
    CS.set_cmap(Cmap)

    CS2 = axes[i].contour(CS, linewidth=6, colors='k', levels=x)

    circleSAN = plt.Circle((100 * cellLength, 100 * cellLength), 70 * cellLength, color='w', fill=False, linewidth=6,
                           linestyle='dashed', alpha=1, zorder=5)
    axes[i].add_artist(circleSAN)

    axes[i].set_title(titles[i])

    if i != len(mapList)- 1:
        axes[i].set_xticklabels([])
        axes[i].set_yticklabels([])
        #axes[i].set_xticks(np.linspace(0,10,20))
        #axes[i].set_yticks(tickLabelsY)
    else:
        axes[i].yaxis.tick_left()
        axes[i].set_xticklabels([])
        axes[i].set_yticklabels([])
        #axes[i].set_xticks(tickLabelsX)
        #axes[i].set_yticks(tickLabelsY)
        #scalebar = ScaleBar(0.01, box_alpha=0.7, pad=1,border_pad=1)
        #axes[i].add_artist(scalebar)


    axes[i].tick_params(length=0, width=0, color='k')
    axes[i].set_aspect('equal')

    numrows, numcols = mapList[i].shape
    def format_coord(x, y):
        col = int(x + 0.5)
        row = int(y + 0.5)
        if col >= 0 and col < numcols and row >= 0 and row < numrows:
            z = mapList[0][row,col]
            return 'x=%1.4f, y=%1.4f, z=%1.4f' % (x, y, z)
        else:
            return 'x=%1.4f, y=%1.4f' % (x, y)

    axes[i].format_coord = format_coord
    #axes[i].clabel(CS2, inline=1, fontsize=10)
    i += 1

fig.subplots_adjust(hspace=0.17, wspace=0.04)


# adding extra white space & placing colorbar there
fig.subplots_adjust(right=0.82, left=0.02, top=0.93)
cbar_ax = fig.add_axes([0.84, 0.232, 0.04, 0.7]) # [0.85, 0.15, 0.05, 0.7]
cbar = fig.colorbar(CS, cax=cbar_ax, ticks=[0, 4, 20, 60, 100, 140])
cbar.set_label('Time, ms')
cbar.set_clim([0., 140.])
cbar.ax.set_yticklabels(np.array([0, 4, 20, 60, 100, 140]))
#cbar.ax.set_yticks()
#plt.tight_layout()
#plt.savefig('chrono_' + str(PROB) + '.eps')
plt.show()
