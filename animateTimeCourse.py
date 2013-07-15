import argparse
#necessary to import os? argparse imports it itself
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def stepData(step, *data):
    i = 0
    print step
    for d in data:
        index = step
        fader = 1.0
        if(index > d.shape[0]-1):
            fader = step/d.shape[0]
            index = d.shape[0]-1
        ax.scatter(d[:index,1], d[:index,2], c=(d[:index,0]/fader), s=8, marker='o', cmap=fadeFromWhiteCmaps[i], vmax=step, lw=0)
        i = i + 1

def animateTimeCourse(inputFiles):
    global ax, fadeFromWhiteCmaps#, data
    fadeFromWhiteCmaps = ['Blues', 'YlOrRd', 'BuGn', 'BuPu', 'YlGn', 'Greens', 'RdPu', 'OrRd', 'Greys', 'Purples']
    data = []
    for f in inputFiles:
        pathToFile = os.path.realpath(f)
        data.append(np.genfromtxt(pathToFile, delimiter=','))
    nSteps = max(d.shape[0] for d in data)
    ylim = max(max(d[:,2]) for d in data)
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    ax.set_xlabel('Fraction voters holding first opinion')
    ax.set_ylabel('Conflicts')
    ax.set_xlim((0,1))
    ax.set_ylim((0,ylim))
    return animation.FuncAnimation(fig, stepData, frames=range(nSteps), fargs=data)
    

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('inputFiles', nargs='+')
    args = parser.parse_args()
    anim = animateTimeCourse(args.inputFiles)
    writer = animation.FFMpegWriter(fps=60, bitrate=8000)
    fileName = args.inputFiles[0]
    us = fileName.find("_")
    fileName = "votingModelAnim" + fileName[us:-4] + ".mkv"
    anim.save(fileName, writer=writer)
