import numpy as np
import os
import matplotlib.pyplot as plt

def plotGraphStats(fileName="voterStats"):
    """Import and plot data from c++ program's
    output"""
    pathToFile = os.path.realpath(fileName)
    graphStats = np.genfromtxt(pathToFile, delimiter = ',')
    plt.figure(1)
    plt.subplot(211)
    plt.plot(graphStats[:,0],graphStats[:,2],'g-')
    plt.subplot(212)
    plt.plot(graphStats[:,1],graphStats[:,2],'b-')
    plt.show()
    
if __name__=='__main__':
    #import sys
    plotGraphStats()
