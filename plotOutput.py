import numpy as np
import os
import matplotlib.pyplot as plt

def plotGraphStats(fileName):
    """Import and plot data from c++ program's
    output"""
    pathToFile = os.path.realpath(fileName)
    toPlot = np.genfromtxt(pathToFile, delimiter = ',')
    if "graphStats" in fileName:
        f, (ax1, ax2) = plt.subplots(2,1)
        ax1.plot(toPlot[:,0],toPlot[:,2],'g-')
        ax1.set_title("1-0 Edges vs. Time")
        #ax1.subplot(212)
        ax2.plot(toPlot[:,1],toPlot[:,2],'b-')
        ax2.set_title("1-0 Edges vs. Opn 1 Vertices")
        plt.show()
    elif "bifData" in fileName:
        toPlot = filter(lambda x: x[3] == 0, toPlot)
        toPlot = np.asarray(toPlot)
        plt.figure(1)
        plt.plot(toPlot[:,0],toPlot[:,1],'og')
        plt.show()
    
if __name__=='__main__':
    import sys
    plotGraphStats(fileName=sys.argv[1])
