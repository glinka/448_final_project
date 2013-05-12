# Script to generate bifData from voterModelAlphaInitDist

from subprocess import call
import os
import numpy as np
import matplotlib.pyplot as plt

for i in range(0,55,5):
    i = i/100.0
    j = 1 - i
    for a in range(0,105,5):
        a = a/100.0
#        sb.Popen([a 2 i j 500], executable=voterModelAlpaInitDist)
#        command = "voterModelAlphaInit " + str(a) + " " + str(2) + " " + str(i) + " " + str(j) + " " + str(500)
        fileLoc = os.getcwd() + "/voterModelAlphaInitDist"
        print fileLoc
        n = 500
        maxIter = 50000
        k = 2
        call([fileLoc, str(a), str(k), str(i), str(j), str(n), str(maxIter)])
bifDataFileName = "bifData_" + str(n) + "_" + str(4) + ".csv"
pathToFile = os.path.realpath(bifDataFileName)
print pathToFile
bifData = np.genfromtxt(pathToFile, delimiter=",")
bifData = filter(lambda x: x[k+1] == 0, bifData)
bifData = np.asarray(bifData)
plt.figure(1)
plt.plot(bifData[:,0],bifData[:,1],'og')
plt.show()
