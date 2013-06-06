import numpy as np
import os
import matplotlib.pyplot as plt

def plotGraphStats(fileName):
    """Import and plot data from c++ program's
    output"""
    pathToFile = os.path.realpath(fileName)
    toPlot = np.genfromtxt(pathToFile, delimiter = ',')
    scaleUp = 16.0/15.0
    scaleDown = 14.0/15.0
    #get graph limits
    l0Low = min(toPlot[:,0])*scaleDown
    l0Up = max(toPlot[:,0])*scaleUp
    l1Low = min(toPlot[:,1])*scaleDown
    l1Up = max(toPlot[:,1])*scaleUp
    l2Low = min(toPlot[:,2])*scaleDown
    l2Up = max(toPlot[:,2])*scaleUp
    nSteps = 10000
    alpha = 0
    plt.hold(True)
    if "graphStats" in fileName:
        f1 = plt.figure(1)
        ax2 = f1.add_subplot(212)
        ax2.plot(toPlot[:,0],toPlot[:,2],'g-')
        ax2.set_title("Number Conflicting Edges vs. Time")
        #ax2.set_ylim([l2Low, l2Up])
	ax2.set_xlim([0,nSteps])
        ax1 = f1.add_subplot(211)
        ax1.plot(toPlot[:,0], toPlot[:,1],'c-')
        ax1.set_title("Minority Fraction vs Time")
	ax1.set_ylim([l1Low,l1Up])
	ax1.set_xlim([0,nSteps])
        f2 = plt.figure(2)
        ax = f2.add_subplot(111)
        #plt.plot(toPlot[:,1],toPlot[:,2],'b-')
        #ax.subplot(212)
        ax.plot(toPlot[:,1],toPlot[:,2],'b-')
        ax.set_title("Fraction Conflicting Edges vs. Fraction Ones")
        ax.set_xlim([l1Low, l1Up])
        ax.set_ylim([l2Low, l2Up])
        plt.show()
    elif "bifData" in fileName:
        toPlot = filter(lambda x: x[2] == 0, toPlot)
        toPlot = np.asarray(toPlot)
        n = toPlot.shape[0]
        #f1 = plt.figure(num=1)
        #ax = f1.add_subplot(111)
        colorList = []
        for i in range(n):
            color = toPlot[i,3]
            if not color in colorList:
                plt.plot(toPlot[i,0], toPlot[i,1],'o',color=[1-2*color, np.cos(np.pi*color), 2*color], label="initial fraction= " + str(color))
                colorList.append(color)
            else:
                plt.plot(toPlot[i,0], toPlot[i,1],'o',color=[1-2*color, np.cos(np.pi*color), 2*color])
        plt.xlabel(r'$\alpha$', fontsize=20)
        plt.ylabel("Final minority fraction")
        plt.legend()
        plt.show()
    
if __name__=='__main__':
    import sys
    plotGraphStats(fileName=sys.argv[1])
