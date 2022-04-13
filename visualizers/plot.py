import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import math as m
import sys
sys.path.insert(0, '/home/lep/projects/results/visualizers')
from common import *


N = 200
map = []
file = sys.argv[1]

with open(file) as f1:
    for line in f1:
        data = line.split()
	for i in range(N**2):        
		map.append(float(data[i]))


lines = [[] for i in range(N)]
map = np.array(map)


left = 0
for j in range(N):	
	right = left + N	
	lines[j] = (map)[left:right]
	left = right


lines = np.vstack(lines)
lines = np.array(lines)

plt.figure()
img = plt.imshow(lines, interpolation = 'none', cmap = 'autumn_r', origin='lower')
#plt.title('Density of LC: PROB = 0.6', fontsize = 25)

fntSize = 20
#plt.xlabel('x, mm', fontsize = fntSize)
#plt.ylabel('y, cells', fontsize = fntSize)
cbar = plt.colorbar()
plt.clim([0,1])
cbar.set_ticks([0, 0.5, 1])
plt.xticks([])
plt.yticks([])
#plt.axis('off')
plt.show()

 

#print lines
