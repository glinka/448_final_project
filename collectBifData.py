# Script to generate bifData from vote

from subprocess import call
import os
import numpy as np

def collectBifData(rewireTo, n, maxIter, collectionInterval):
	for i in range(5,55,5):
	    i = i/100.0
	    j = 1 - i
	    for a in range(0,120,20):
		a = a/100.0
		fileLoc = os.getcwd() + "/vote"
		k = 2
		call([fileLoc, rewireTo, str(n), str(maxIter), str(collectionInterval), str(a), str(k), str(i), str(j)])

if __name__=="__main__":
	import sys
	collectBifData(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])






