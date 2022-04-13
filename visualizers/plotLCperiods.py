#!/usr/bin/env pytnon
# vim: set fileencoding=utf-8 ts=4 sw=4 expandtab:

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib
import sys
from common import *

'''
# russian letters
from matplotlib import rc
rc('font',**{'family':'sansserif'})
rc('text', usetex=True)
rc('text.latex',unicode=True)
rc('text.latex',preamble='\usepackage[utf8]{inputenc}')
rc('text.latex',preamble='\usepackage[russian]{babel}')
'''

font = {'family': 'normal', 'weight': 'normal', 'size': 40}
matplotlib.rc('font', **font)

cellLength = 0.07

def ReadCSV(fileName):
	# reading 1st line
	f = open(fileName)
	names = f.readline()
	data = np.loadtxt(fileName, delimiter = ',', skiprows = 1, unpack = True)
	return names, data


def Concatentate(t, timeArrayInBeginning, timeArrayWhole, dataInBeginning, dataWhole):
	
	timeArrayInBeginning = list(timeArrayInBeginning)
	timeArrayWhole = list(timeArrayWhole)

	dataInBeginning = list(dataInBeginning)
	dataWhole = list(dataWhole)

	dataAfter_t = []
	timeArrayAfter_t = []

	for i in range(len(timeArrayWhole)):
		#print timeArrayWhole[i]
		if timeArrayWhole[i] > t:
			dataAfter_t.append(dataWhole[i])
			timeArrayAfter_t.append(timeArrayWhole[i])

	dataWholeConcatanated = dataInBeginning + dataAfter_t	
	timeArrayConcatanated = timeArrayInBeginning + timeArrayAfter_t
	
	print timeArrayAfter_t
	
	return np.array(timeArrayConcatanated), np.array(dataWholeConcatanated)




# for concatenation
files1stSec = ['/media/lep/Elements/Chrono_copy_big_files/PROB00/E_before_ACH_0_20_secs.dat_less_1sec.txt', '/media/lep/Elements/Chrono_copy_big_files/PROB04/potentials_PROB=0.4.dat_less_1sec.txt', \
				'/media/lep/Elements/Chrono_copy_big_files/PROB06/potentials_PROB=0.6.dat_less_1sec.txt']


numberOfFiles = int(sys.argv[1])

files = []

#reading files
for i in range(2, numberOfFiles + 2):
	files.append(sys.argv[i])


t, numberOfLC, periods, smth = [], [], [], []
for i in range(numberOfFiles):
	names, dataWhole = ReadCSV(files[i])
	names2, dataInBeginning = ReadCSV(files1stSec[i])

	timeConcatenated, dataConcatenated1 = Concatentate(1., dataInBeginning[0], dataWhole[0], dataInBeginning[1], dataWhole[1])
	timeConcatenated, dataConcatenated2 = Concatentate(1., dataInBeginning[0], dataWhole[0], dataInBeginning[2], dataWhole[2])

	#timeConcatenated, dataConcatenated3 = Concatentate(1., dataInBeginning[0], dataWhole[0], dataInBeginning[3], dataWhole[3])

	t.append(timeConcatenated)
	numberOfLC.append(dataConcatenated1)
	periods.append(dataConcatenated2)
	#smth.append(data[3])

names = names.split('||')
colors = ['r-', 'g-', 'b-']
labels = ['P = 0.0', 'P = 0.4', 'P = 0.6']
#labels = ['500 frames', '10 frames']

#plotting
plt.figure()
for i in range(numberOfFiles):
	plt.plot(t[i], periods[i], colors[i], linewidth = 6, label = labels[i], markersize = 10, markeredgecolor='none')




leg = plt.legend(loc='upper right', prop={'size':40}, markerscale = 3, numpoints=1)
for l in leg.legendHandles:            
	l.set_linewidth(10)


#plt.gca().set_ylim([0.0, 65])
#plt.legend(numpoints=1)
plt.axvline(x = 60, linestyle= '--', linewidth=6, color='k')
plt.xlim([0., 1.])
plt.ylabel('Cycle length, s')
plt.xlabel('Time, s')
plt.xticks([30, 60, 90, 120])
#plt.yticks([0, 30, 60, 90, 120])
plt.grid('on')
plt.show()

