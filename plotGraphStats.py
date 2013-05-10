import numpy as np
import os
import matplotlib.pyplot as plt

def plotGraphStats(fileName):
    """Import and plot data from c++ program's
    output"""
    pathToFile = os.path.realpath(fileName)
    graphStats = np.genfromtxt(pathToFile, delimiter = ',')
    f, (ax1, ax2) = plt.subplots(2,1)
    ax1.plot(graphStats[:,0],graphStats[:,2],'g-')
    ax1.set_title("1-0 Edges vs. Time")
    #ax1.subplot(212)
    ax2.plot(graphStats[:,1],graphStats[:,2],'b-')
    ax2.set_title("1-0 Edges vs. Opn 1 Vertices")
    plt.show()
    
if __name__=='__main__':
    import sys
    plotGraphStats(fileName=sys.argv[1])
