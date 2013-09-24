import numpy as np

def get_conflicts(adj_matrix, opns_vector):
    n = opns_vector.shape[0]
    conflicts = 0
    for i in range(n):
        for j in range(i+1, n):
            conflicts = conflicts + adj_matrix[i,j]*np.absolute(opns_vector[i] - opns_vector[j])
    return conflicts

def get_minority_fraction(adj_matrix, opns_vector):
    n = opns_vector.shape[0]
    unique_opns = {}
    for i in range(n):
        if str(opns_vector[i]) in unique_opns.keys():
            unique_opns[str(opns_vector[i])] = unique_opns[str(opns_vector[i])] + 1
        else:
            unique_opns[str(opns_vector[i])] = 1
    min = n
    for val in unique_opns.itervalues():
        if val < min:
            min = val
    return 1.0*min/n
