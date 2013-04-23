# Voting model implementation (based on Shi, Mucha and Durrett's 
# "A multi-opinion evolving voter model with infinitley many phase transistions")

import numpy as np
import matplotlib.pyplot as plt

def vote(n=1000, avgDeg=4.0, u=(0.4,0.6), a=0.8, rewireTo="random", maxIter=10000, timeInterval=10):
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
    #print A, Opns, conflicts
    iters = 0
    N10timeCourse = np.zeros(int(maxIter/timeInterval))
    N1timeCourse = np.zeros(int(maxIter/timeInterval))
    stepTimeCourse = np.zeros(int(maxIter/timeInterval))
    step = 0
    if rewireTo == 'random':
        while (conflicts > 0) & (iters < maxIter):
            #!!! could also choose edge, may reduce compuation time !!!
            chosenVertex = int(np.floor(n*np.random.random_sample(1)))
            #discard vertices until deg(v) != 0
            #could potentially improve in two ways:
            #1. keep track of vertex degrees in separate array (easy)
            #2. somehow pop out vertices that have degree zero, or at least
            #remove them as they are found in the below iteration (?)
            while sum(A[chosenVertex][:]) == 0:
                chosenVertex = int(np.floor(n*np.random.random_sample(1)))
            #generate list of adjacent vertices, V
            V = []
            for j in range(n):
                if A[chosenVertex][j] != 0:
                    V.append(j)
            #V = np.array(V)
            numberAdj = len(V) #V.size
            #!!! Previous implementations allowed the relocation of the edge
            #even when Opns[v1] = Opns[v2]; however, this resulted in significantly
            #different system dynamics than those described in the paper, in which
            #no actions of any kind are taken if Opns[v1] = Opns[v2] !!!
            neighbor = V.pop(int(np.floor(numberAdj*np.random.random_sample(1))))
            numberAdj = numberAdj - 1
            while (numberAdj > 0) & (Opns[chosenVertex] == Opns[neighbor]):
                neighbor = V.pop(int(np.floor(numberAdj*np.random.random_sample(1))))
                numberAdj = numberAdj - 1
            if Opns[chosenVertex] != Opns[neighbor]:
                actionToPerform = np.random.random_sample(1)
                conflictCounter = 0
                if actionToPerform > a:
                    #force chosenVertex to agree with neighbor
                    Opns[chosenVertex] = Opns[neighbor]
                    #update conflicts
                    for j in range(n):
                        if (A[chosenVertex][j] != 0) & (Opns[j] != Opns[chosenVertex]):
                            conflictCounter = conflictCounter + 1
                        elif (A[chosenVertex][j] != 0):
                            conflictCounter = conflictCounter - 1
                    conflicts = conflicts + conflictCounter
                else:
                    #log changes in conflicts, remove edge and add new edge,
                    #non-parallel, non-loop, to chosenVertex
                    if Opns[chosenVertex] != Opns[neighbor]:
                        conflictCounter = -1
                    A[chosenVertex][neighbor] = 0
                    A[neighbor][chosenVertex] = 0
                    newNeighbor = int(np.floor(n*np.random.random_sample(1)))
                    #check the added edge will not be a loop or in parallel with a previously existing one
                    #this loop will require many iterations to exit if deg(v) ~ O(n), potential slowdown
                    while (A[chosenVertex][newNeighbor] != 0) | (newNeighbor == chosenVertex):
                        newNeighbor = int(np.floor(n*np.random.random_sample(1)))
                    A[chosenVertex][newNeighbor] = 1
                    A[newNeighbor][chosenVertex] = 1
                    if Opns[chosenVertex] != Opns[newNeighbor]:
                        conflictCounter = conflictCounter + 1
                    conflicts = conflicts + conflictCounter
            if iters % timeInterval == 0:
                step = iters/timeInterval
                graphStats = calcGraphStatistics(A, Opns, len(u))
                N1timeCourse[step] = graphStats[1][0]
                N10timeCourse[step] = graphStats[2]
                stepTimeCourse[step] = step + 1
            iters = iters + 1
            #print conflicts == calcConflict(A, Opns)
            #plot results
    plt.figure(1)
    plt.subplot(211)
    plt.plot(N1timeCourse[:step],N10timeCourse[:step],'g-')
    plt.subplot(212)
    plt.plot(stepTimeCourse[:step],N10timeCourse[:step])
    plt.show()
    return A, Opns, p

#!!! consider making truly general for n opinions
def calcGraphStatistics(A, Opns, numOpns):
    totalEdges = 0
    opnFractions = np.zeros(numOpns)
    discordantEdges = 0
    n = A.shape[0]
    for i in range(n):
        currentOpn = int(Opns[i])
        opnFractions[currentOpn - 1] = opnFractions[currentOpn - 1] + 1
        for j in range(i+1,n):
            totalEdges = totalEdges + A[i][j]
            if (A[i][j] != 0) & (Opns[j] != currentOpn):
                discordantEdges = discordantEdges + 1
    opnFractions = [(1.0*opnFractions[i])/n for i in range(numOpns)]
    discordantEdges = discordantEdges*1.0/totalEdges
    return (totalEdges, opnFractions, discordantEdges)

def calcConflict(A, Opns):
    conflicts = 0
    n = A.shape[0]
    for i in range(n):
        currentOpn = Opns[i]
        for j in range(i+1,n):
            if (A[i][j] != 0) & (Opns[j] != currentOpn):
                conflicts = conflicts + 1
    return conflicts



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


    
    
    
