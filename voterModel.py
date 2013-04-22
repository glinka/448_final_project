# Voting model implementation (based on Shi, Mucha and Durrett's 
# "A multi-opinion evolving voter model with infinitley many phase transistions")

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import lambertw

def vote(n=10, avgDeg=4.0, u=(0.5,0.5), a=0.5, rewireTo="random"):
    """simulate voting model with k opinions, alpha = a
       and rewiring scheme of rewireTo (with default parameters
       2, 0.5 and 'random', respectively)""" 
    #p = probability of edge in grap
    #v = current vertex in iteration (v in [1,n])
    #w = number of "steps" till next edge, chosen from geometric
    #distribution as in Physical Review E 71, 036113 (2005)
    #A is adjacency matrix
    #O is corresponding opinion matrix
    #!! Combine A and O into matrix of dimension (2,n,n)? !!
    if sum(u) != 1:
        return 'invalid probability distribution, sum(u) != 1'
    p = avgDeg/(n-1)
    v = 1
    w = -1
    A = np.zeros((n,n))
    Opns = np.zeros((n,1))
    while v <= n:
        #r determines where the next vertex is added
        #o determines the opinion of that vertex
        #w calculates the number of vertices to skip
        #based on r and p
        rv1 = np.random.random_sample(1)
        w = w + 1 + int(np.floor(np.log(1-rv1)/np.log(1-p)))
        while (w >= v) & (v <= n):
            rv2 = np.random.random_sample(1)
            partialSum = 0.0
            indexCounter = 0
            while partialSum < rv2:
                partialSum = partialSum + u[indexCounter]
                indexCounter = indexCounter + 1
            Opns[v-1][0] = indexCounter
            w = w - v
            v = v + 1
        if v < n:
            A[v-1][w-1] = 1
            A[w-1][v-1] = 1
    #calculate number of disagreeing vertices
    conflicts = 0
    for i in range(n):
        currentOpinion = Opns[i]
        for j in range(i+1,n):
            if currentOpinion + Opns[j] != 2*currentOpinion:
                conflicts = conflicts + 1
    print conflicts
    return A, Opns, p
            






def checkDegrees(A, p):
    """ensures A is properly initialized, in regards to average
    degree of each vertex"""
    dim = np.shape(A)[0]
    degrees = np.zeros(dim)
    avg = 0.0
    for i in range(dim):
        for j in range(dim):
            degrees[i] = degrees[i] + A[i][j]
        avg = avg + degrees[i]
    avg = avg / dim
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    n, bins, patches = ax1.hist(degrees, bins=dim/2)
    model = np.random.binomial(dim, p, np.ceil(dim**2*p/2))
    ax2 = fig.add_subplot(212)
    n, bins, patches = ax2.hist(model, bins=dim/2)
    plt.show()
    return avg, degrees


    
    
    
