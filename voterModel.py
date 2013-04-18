# Voting model implementation (based on Shi, Mucha and Durrett's 
# "A multi-opinion evolving voter model with infinitley many phase transistions")

import numpy as np

def vote(k = 2, a = 0.5, rewireTo = "random", ):
    """simulate voting model with k opinions, alpha = a
       and rewiring scheme of rewireTo (with default parameters
       2, 0.5 and 'random', respectively)""" 
    
