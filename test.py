import numpy as np
import matplotlib.pyplot as plt

def test():
    nfns = 2
    ndata = 50
    var = 0.5
    A = np.zeros([ndata,2])
    x = np.linspace(1, 10, ndata) + var*np.random.randn(ndata)
    y = -1*np.linspace(11 ,20, ndata) + var*np.random.randn(ndata)
    for i in range(ndata):
        A[i, 0] = x[i]
        A[i, 1] = 1
    coeffsA = np.dot(np.linalg.inv(np.dot(np.transpose(A), A)), np.dot(np.transpose(A), y))
    coeffsB = np.dot(np.dot(np.linalg.inv(np.dot(np.transpose(A), A)), np.transpose(A)), y)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, coeffsA[0]*x+coeffsA[1])
    ax.plot(x, coeffsB[0]*x+coeffsB[1])
    ax.plot(x, y)
    plt.show()

if __name__=="__main__":
    test()
