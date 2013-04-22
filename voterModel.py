# Voting model implementation (based on Shi, Mucha and Durrett's 
# "A multi-opinion evolving voter model with infinitley many phase transistions")

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import lambertw

def vote(n=10, h=4.0, k=2, a=0.5, rewireTo="random"):
    """simulate voting model with k opinions, alpha = a
       and rewiring scheme of rewireTo (with default parameters
       2, 0.5 and 'random', respectively)""" 
    #create random vector uniformly distributed between [0,1) of size n
    v = np.array(np.random.random_sample(n), ndmin=2)
    vT = np.array(np.random.random_sample(n), ndmin=2).T
#    vT = v.T
    A = np.dot(vT,v)
    #probability of edge in graph, threshold A with p
    p = h/(n-1)
#    print p, v, A
    thresholdConst = np.real(np.exp(lambertw((p-1)/np.exp(1),k=-1)+1))
    threshold = np.vectorize(lambda x: np.ceil(x-thresholdConst))
    A = threshold(A)
    for i in range(np.shape(A)[0]):
        A[i][i] = 0;
    return A, p

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


    
    
    
