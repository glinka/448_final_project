# Script to generate bifData from voterModelAlphaInitDist

from subprocess import call
import os
import numpy as np

def collectBifData(n, maxIter, rewireTo="random"):
	for i in range(0,55,5):
	    i = i/100.0
	    j = 1 - i
	    for a in range(0,110,20):
		a = a/100.0
		fileLoc = os.getcwd() + "/voterModelAlphaInitDist"
		k = 2
		call([fileLoc, str(a), str(k), str(i), str(j), str(n), str(maxIter), rewireTo])

if __name__=="__main__":
	import sys
	collectBifData(sys.argv[1], sys.argv[2], sys.argv[3])

