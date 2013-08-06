import numpy as np
import os
import matplotlib.pyplot as plt

def plotGraphStats(fileNames):
    """Import and plot data from c++ program's
    output"""
    nDataFig = plt.figure(len(fileNames)+1)
    nDataAx = nDataFig.add_subplot(111)
    nData = 0
    i = 0
    for fileName in fileNames:
        pathToFile = os.path.realpath(fileName)
        toPlot = np.genfromtxt(pathToFile, delimiter = ',')
        if "graphStats" in fileName:
            foundUnderScores = 0
            alpha = 0
            while foundUnderScores < 3:
                alpha = fileName.find("_", alpha+1)
                foundUnderScores = foundUnderScores + 1
            alpha2 = fileName.find("_", alpha+1)
            alpha = float(fileName[alpha+1:alpha2])
            nData = nData + 1
            i = i + 1
            scaleUp = 16.0/15.0
            scaleDown = 14.0/15.0
            #get graph limits
            l0Low = min(toPlot[:,0])*scaleDown
            l0Up = max(toPlot[:,0])*scaleUp
            l1Low = min(toPlot[:,1])*scaleDown
            l1Up = max(toPlot[:,1])*scaleUp
            l2Low = min(toPlot[:,2])*scaleDown
            l2Up = max(toPlot[:,2])*scaleUp
            nSteps = toPlot.shape[0]
            f1 = plt.figure(i)
            ax2 = f1.add_subplot(211, ylabel="Number conflicting edges")
            ax2.plot(toPlot[:,0],toPlot[:,2], '-', color="#47D147")
            ax2.set_ylim([l2Low, l2Up])
            ax2.set_xlim([0,nSteps])
            ax1 = f1.add_subplot(212, sharex=ax2, ylabel="Minority fraction", xlabel="Time")
            ax1.plot(toPlot[:,0], toPlot[:,1], '-', color="#70AAFF")
            ax1.set_ylim([l1Low,l1Up])
            plt.setp(ax2.get_xticklabels(), visible=False)
            plt.subplots_adjust(hspace=0)
            saveFilename = fileName[:-4] + "_fig1" + ".png"
            plt.savefig(saveFilename, edgecolor='k')
            nDataAx.plot(toPlot[:,1],toPlot[:,2],'b-', label=r'$\alpha$' + " = " + str(alpha),  color=[1-alpha, np.cos(np.pi*alpha/2), alpha])
        elif "bifData" in fileName:
            i = i + 1
            fig = plt.figure(i)
            ax = fig.add_subplot(111)
            alpha = 0
            #filtering by completion of runs: x[2] == 0
            toPlot = filter(lambda x: ((x[2] == 0) & (x[3] != 0)), toPlot)
            toPlot = np.asarray(toPlot)
            n = toPlot.shape[0]
            #ids is a 2d list, with the first holding
            #the different values of alpha, the second
            #the corresponding initial minority fraction
            ids = {}
            colorList = []
            nu0 = 0
            nAlpha = 0
            for i in range(n):
                alpha = int(100*toPlot[i,0])
                u0 = int(100*toPlot[i,3])
                uf = toPlot[i,1]
                if not alpha in ids:
                    ids[alpha] = {}
                if not u0 in ids[alpha]:
                    ids[alpha][u0] = []
                ids[alpha][u0].append(uf)
                nAlpha = len(ids)
                if len(ids[alpha]) > nu0:
                    nu0 = len(ids[alpha])
            toPlot = np.zeros((nu0,nAlpha,3))
            u0s = []
            u0f = []
            j = 0
            for a in ids:
                u0s = ids[a].keys()
                u0s.sort()
                i = 0
                for u0 in u0s:
                    ids[a][u0] = np.average(np.asarray(ids[a][u0]))
                    toPlot[i,j,0] = a/100.0
                    toPlot[i,j,1] = ids[a][u0]
                    i = i + 1
                if len(u0s) > len(u0f):
                    u0f = u0s
                j = j + 1
            for i in range(nu0):
                u0 = u0f[i]/100.0
                ax.plot(toPlot[i,:,0], toPlot[i,:,1], 'o', color=[1-2*u0, np.cos(np.pi*u0), 2*u0], label="initial fraction= " + str(u0))
            ax.set_xlabel(r'$\alpha$', fontsize=20)
            ax.set_ylabel("Final minority fraction")
            colorList.sort()
            ax.legend(loc=3)
            ax.set_xlim([0,1])
            ax.set_ylim([0,0.5])
            saveFilename = fileName[:-4] + ".png"
            plt.savefig(saveFilename, edgecolor='k')
    nDataAx.set_title("Number of conflicting edges vs. Minority fraction")
    nDataAx.legend(loc=2)
    saveFilename = "nData.png"
    plt.savefig(saveFilename, edgecolor='k')
    plt.show()
    
if __name__=='__main__':
    import sys
    fileNames = []
    plotGraphStats(sys.argv[1:])
