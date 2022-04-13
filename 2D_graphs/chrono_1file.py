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
#import seaborn as sbn


import matplotlib.dates as mdates
import matplotlib
import sys
import os
import numba as nb
from matplotlib.ticker import FormatStrFormatter

from matplotlib import rcParams
#rcParams.update({'figure.autolayout': True})

font = {'family': 'normal', 'weight': 'normal', 'size': 40}
matplotlib.rc('font', **font)
import math as m


dt= 1e-1
SZ = 40

'''
Cmap = []
for i in range(256):
    Cmap.append(cm.autumn(i))

Cmap = np.array(Cmap)



Cmap = TransformColormap(rpd, Cmap)
Cmap = crs.LinearSegmentedColormap.from_list('my_colormap', Cmap, 256)
print Cmap
'''



dirName = '/media/lep/Elements/Chrono_copy2/PROB00/'
files = ['f11.dat', 'f12.dat', 'f13.dat', 'f21.dat', 'f22.dat', 'f23.dat']
#files = ['v11.dat', 'v12.dat', 'v13.dat', 'v21.dat', 'v22.dat', 'v23.dat']
titles = ['a', 'b', 'c', 'd', 'e', 'f']


fileList = [dirName + item for item in files]

print fileList
mapList = []

for file in fileList:
    mapList.append(np.reshape(np.fromfile(file, dtype=np.float32, count=-1, sep=""), (200,200)))

# processing nans
'''
for k in range(200):
    for m in range(200):
        if np.isnan(mapList[3][k, m]):
            mapList[3][k, m] = 0.
'''


print mapList[3]

Levels=[np.arange(1492,1527,5), np.arange(1527,1577,10), np.arange(1597,1697,20), np.arange(1737,3000,40)]
Levels=list(it.chain.from_iterable(Levels))
Levels= (np.array(Levels)- 1492)*dt
print 'Levels = ', Levels


cellLength =0.07
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


fig, axes = plt.subplots(1, 1)
fig.subplots_adjust(hspace = 0.2, wspace=0.2)
#axes = axes.ravel()

Cmap = []
i = 0

mapList[i][:,:] = (mapList[i][:,:] - np.amin(np.amin(mapList[i])))*dt
print mapList[i]+

#mapList[i] = gaussian_filter(mapList[i], sigma=0) # gaussian filter looses critical details

CS = axes.contourf(mapList[i], levels=Levels, extent=[0, 14, 0,14])
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

    CS2 = axes.contour(CS, linewidth=6, colors='k', levels = x)
    circleSAN = plt.Circle((100 * cellLength, 100 * cellLength), 70 * cellLength, color='w', fill=False, linewidth=6,
                           linestyle='dashed', alpha=1, zorder=5)
    axes.add_artist(circleSAN)

    axes.set_title(titles[i], y=1.03)

    numrows, numcols = mapList[i].shape
    def format_coord(x, y):
        col = int(x + 0.5)
        row = int(y + 0.5)
        if col >= 0 and col < numcols and row >= 0 and row < numrows:
            z = mapList[0][row,col]
            return 'x=%1.4f, y=%1.4f, z=%1.4f' % (x, y, z)
        else:
            return 'x=%1.4f, y=%1.4f' % (x, y)

    axes.format_coord = format_coord
    #axes[i].clabel(CS2, inline=1, fontsize=10)
    i += 1

# adding extra white space & placing colorbar there
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cbar = fig.colorbar(CS, cax=cbar_ax, ticks=[0, 4, 20,  60, 100, 140])
cbar.set_label('ms')
#cbar.set_clim([0, 150])
#cbar.set_ticks(np.linspace(0, 150, 10))
#cbar.ax.set_yticks([0, 1400])

plt.show()
