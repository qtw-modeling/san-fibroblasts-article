from common import *

cellXcoord = 100
cellYcoord = 100
N = 200
numCellsTotal = N**2
outputStep = 5.0e-3


def TimeFromFrame(frame):
    return frame/outputStep


def FrameFromTime(time):
    return int(time / outputStep) 


def CalculateAmplitude(array1D, timeStart, timeEnd):
    array = array1D[0:FrameFromTime(timeEnd - timeStart)]  # timeACH --- time of superfusion of ACh
    array = np.absolute(array)
    return array.max()


# reading files
numFiles = int(sys.argv[1])
fileName = []

for i in range(2, numFiles + 2):
    fileName.append(sys.argv[i])


def Idx(i, j, frame=0):
    return frame * numCellsTotal + i + N * j


potential = [[] for i in range(numFiles)]
# pacemaker cell at (107, 93)
cellIdxList = [Idx(100, 100), Idx(100, 170), Idx(100, 190)]

initialTime = 0
endTime = 2

#print FrameFromTime(ACHtime)

[initialFrame, endFrame] = TimeFromFrame(np.array([initialTime, endTime]))

#numOutputFrames = 
#endFrame = initialFrame + numOutputFrames
#step = int((endFrame - initialFrame) / numOutputPoints)
outputFrames = range(int(initialFrame), int(endFrame + 1))
sizeofFormat = 4 # in bytes 

#print frames
ytikz = np.linspace(-90, 40, 3)

for i in range(numFiles):
	f = open(fileName[i], 'rb')
	
    	for frame in outputFrames:
	 	f.seek(frame * numCellsTotal * sizeofFormat, os.SEEK_SET)
	    	potentialsWholeMap = np.fromfile(f, dtype = np.dtype('f'), count = numCellsTotal)
	    	potential[i].append(potentialsWholeMap[cellIdxList[2]])
	
    	potential[i] = np.array(potential[i])

    	frames = np.array(outputFrames)
    	time = np.array(frames * outputStep)
    	fig, ax = plt.subplots()
    	plt.plot(time, potential[i], 'k', linewidth = 4)

    	#plt.axvline(x = 60, linestyle= '--', linewidth=6, color='k')
    	plt.xlabel('Time, s')
    	plt.ylabel('[], moles/l')
        #plt.yticks(ytikz)
        #plt.ylim([ytikz[0], ytikz[2]])
	#ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
	#plt.xticks(np.linspace(60, 64, 3))
	plt.grid('on')
    	forOneFile = sys.argv[2] + '.png'
    	#plt.savefig(forOneFile, format = 'png', dpi = 200)
	f.close()
plt.show()

