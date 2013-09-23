import numpy as np
import os
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid.inset_locator as insetAxes

def plotGraphStats(fileNames):
    """Import and plot data from c++ program's
    output"""
    nDataFig = plt.figure(len(fileNames)+1)
    nDataAx = nDataFig.add_subplot(111)
    nData = 0
    haveGraphStats = False
    plotIndividuals = True
    i = 0
    for fileName in fileNames:
        pathToFile = os.path.realpath(fileName)
        toPlot = np.genfromtxt(pathToFile, delimiter = ',')
        if "convergenceData" in fileName:
            #runs per step
            rps = 10
            nSteps = toPlot.shape[0]/rps
            fig = plt.figure();
            ax = fig.add_subplot(111);
            for i in range(nSteps):
                if toPlot[i*rps, 1] == 0.482:
                    toPlot[i*rps:(i+1)*rps, 1] = 0
            xData = [sum(toPlot[i*rps:(i+1)*rps, 1])/rps for i in range(nSteps)]
            yData = [sum(toPlot[i*rps:(i+1)*rps, 0])/rps for i in range(nSteps)]
            xMax = 1.0*np.max(xData)
            yMax = 1.0*np.max(yData)
            xUpLim = xMax
            yUpLim = 1000000
            ax.scatter(xData, yData, c='r', s=36)
            ax.set_xlabel('Projection interval')
            ax.set_ylabel('Steps to consensus (log scale)')
            ax.set_yscale('log')
            ax.set_ylim((10000, yUpLim))
            ax.set_xlim((0, xMax))
            insets = []
            #copy list to allow modifications
            fileList = fileNames
            fileList.remove(fileName)
            for i in range(nSteps):
                for f in fileList:
                    #if have graphstats with same projection step, plot in subplot
                    if "graphStats" and str(xData[i]) in f:
                        pathToFile = os.path.realpath(f)
                        toPlotInset = np.genfromtxt(pathToFile, delimiter = ',')
                        insets.append(insetAxes.inset_axes(ax, width="50%", height="50%", loc=3, bbox_to_anchor=(xData[i]/xMax, np.log10(yData[i])/np.log10(yUpLim), 0.4, 0.4), bbox_transform=ax.transAxes))
                        insets[i].set_xticks([])
                        insets[i].set_yticks([])
                        insets[i].plot(toPlotInset[:,1], toPlotInset[:,2])
                        break
            plt.show(fig)
        #don't do individual plotting right now
        if ("graphStats" in fileName) & (plotIndividuals):
            haveGraphStats = True
            paramstr = fileName
            us = paramstr.find("_")
            paramstr = paramstr[us+1:]
            params = {}
            while us > 0:
                us = paramstr.find("_")
                us2 = paramstr.find("_", us+1)
                params[paramstr[:us]] = float(paramstr[us+1:us2])
                if "initDist" in paramstr[:us]:
                    period = paramstr.find(".", us2+1)
                    period2 = paramstr.find(".", period+1)
                    params[paramstr[:us]] = [float(paramstr[us+1:us2]), float(paramstr[us2+1:period2])]
                    us = -1
                else:
                    paramstr = paramstr[us2+1:]
            for key in params:
                if("initDist" not in key):
                    if(params[key] % 1 == 0):
                        params[key] = int(params[key])
            print params
            alpha = params['alpha']
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
    if haveGraphStats:
        nDataAx.set_title("Number of conflicting edges vs. Minority fraction")
        nDataAx.legend(loc=2)
        saveFilename = "nData.png"
        plt.savefig(saveFilename, edgecolor='k')
        plt.show()
    else:
        plt.close(nDataFig)
    
if __name__=='__main__':
    import sys
    fileNames = []
    plotGraphStats(sys.argv[1:])
