# Voting model implementation (based on Shi, Mucha and Durrett's 
# "A multi-opinion evolving voter model with infinitley many phase transistions")

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import lambertw

def vote(n=10, avgDeg=4.0, u=(0.5,0.5), a=0.5, rewireTo="random"):
    """simulate voting model with k opinions, alpha = a
       and rewiring scheme of rewireTo (with default parameters
       2, 0.5 and 'random', respectively)""" 
    #p = probability of edge in graph
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
        #rv1 determines where the next edge is added
        #rv2 determines opinion of each vertex, based on 
        #input probabilities
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
            if (A[i][j] != 0) & (Opns[j] != currentOpinion):
                conflicts = conflicts + 1
    print A, Opns, conflicts
    maxIter = 100000
    iters = 0
    if rewireTo == 'random':
        while (conflicts > 0) & (iters < maxIter):
            #!! could also choose edge, may reduce compuation time !!
            #!! implement fn: conflictCalc(A)?? !!
            chosenVertex = int(np.floor(n*np.random.random_sample(1)))
            actionToPerform = np.random.random_sample(1)
            #discard vertices until deg(v) != 0
            while sum(A[chosenVertex][:]) == 0:
                chosenVertex = int(np.floor(n*np.random.random_sample(1)))
            #generate list of adjacent vertices, V
            V = []
            for j in range(n):
                if A[chosenVertex][j] != 0:
                    V.append(j)
            V = np.array(V)
            numberAdj = V.size
            neighbor = V[int(np.floor(numberAdj*np.random.random_sample(1)))]
            conflictCounter = 0
            if actionToPerform > a:
                if Opns[chosenVertex] != Opns[neighbor]:
                    Opns[chosenVertex] = Opns[neighbor]
                    for j in range(n):
                        if (A[chosenVertex][j] != 0) & (Opns[j] != Opns[chosenVertex]):
                            conflictCounter = conflictCounter + 1
                        elif (A[chosenVertex][j] != 0):
                            conflictCounter = conflictCounter - 1
                conflicts = conflicts + conflictCounter
            else:
                if Opns[chosenVertex] != Opns[neighbor]:
                    conflictCounter = -1
                A[chosenVertex][neighbor] = 0
                A[neighbor][chosenVertex] = 0
                neighbor = int(np.floor(n*np.random.random_sample(1)))
                #check the added edge will not be a loop or in parallel with a previously existing one
                while (A[chosenVertex][neighbor] != 0) | (neighbor == chosenVertex):
                    neighbor = int(np.floor(n*np.random.random_sample(1)))
                A[chosenVertex][neighbor] = 1
                A[neighbor][chosenVertex] = 1
                if Opns[chosenVertex] != Opns[neighbor]:
                    conflictCounter = conflictCounter + 1
                conflicts = conflicts + conflictCounter
            iters = iters + 1
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


    
    
    
