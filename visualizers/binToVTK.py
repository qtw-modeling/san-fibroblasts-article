import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import sys
import os


fileName = sys.argv[1]




nx, ny, nz = 200, 200, 1 
numCellsTotal = nx * ny * nz
outputStep = 5.0e-3




def To2darray(map, numCells):
	lines = [[] for i in range(numCells)]
	map = np.array(map)

	left = 0
	for j in range(numCells):	
		right = left + numCells
		lines[j] = (map)[left:right]
		left = right

	lines = np.vstack(lines)
	lines = np.array(lines)
	return lines



def writeVTK(array, frame):
	fOut = open('ap_%d.vtk' %  frame, 'w+')
	# VTK output
    	fOut.write('# vtk DataFile Version 3.0\nSolution\nASCII\nDATASET RECTILINEAR_GRID\nDIMENSIONS %d %d %d\n' % (nx + 1, ny + 1, nz))

    
    	fOut.write('X_COORDINATES %d double\n' % (nx + 1))

    	for i in range(nx + 1):
        	fOut.write('%d ' % i)
        
    	fOut.write('\n')
   
    	fOut.write('Y_COORDINATES %d double\n' % (ny + 1))
    	for i in range(ny + 1):
    		fOut.write('%d ' % i)
       
    	fOut.write('\n')

    	fOut.write('Z_COORDINATES 1 double\n')
	fOut.write('0')
	
	fOut.write('\n') 

    	fOut.write('CELL_DATA %d\n' % (nx * ny))   
    	# printing all the data in the file
    	fOut.write('SCALARS type double\nLOOKUP_TABLE default\n')
      
    	for idx in range(nx * ny):
    		fOut.write('%e ' % array[idx])
          
        	if ((idx + 1)%nx == 0):
           		fOut.write('\n')

    	fOut.close()




f = open(fileName, 'rb')



def TimeFromFrame(frame):
    return frame/outputStep


def FrameFromTime(time):
    return int(time / outputStep) 


initialTime = 20
endTime = 70

#print FrameFromTime(ACHtime)

[initialFrame, endFrame] = TimeFromFrame(np.array([initialTime, endTime]))

#numOutputFrames = 
#endFrame = initialFrame + numOutputFrames
#step = int((endFrame - initialFrame) / numOutputPoints)
outputFrames = range(int(initialFrame), int(endFrame + 1))
sizeofFormat = 4 # in bytes 


for frame in outputFrames:
    f.seek(frame * numCellsTotal * sizeofFormat, os.SEEK_SET)
    potentialsWholeMap = np.fromfile(f, dtype = np.dtype('f'), count = numCellsTotal)
    
    
    #out2Darray = To2darray(potentialsWholeMap, N)
    writeVTK(potentialsWholeMap, frame)

    if (frame % 100 == 0):
		print '#frame', frame
	
print 'Finished making VTKs from Binary'
 


