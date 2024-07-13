'''
Utility functions
'''

import csv
import numpy as np


def newtonRaphson(y, dydx, init_guess, tol, maxiter=1000):
    '''
    Calculate root of single variable function using Newton-Raphson method
    '''
    
    x = init_guess

    for i in range(maxiter):
        x_new = x - y(x) / dydx(x)
        if np.abs(x_new - x) < tol:
            return x_new
        x = x_new
    raise RuntimeError("newtonRaphson: This function did not converge to a solution.")
    
def norm(vec):
    '''
    Vector norm
    '''
    return np.linalg.norm(vec)

def acos(x):
    '''
    Inverse Cosine
    '''
    return np.arccos(np.clip(x, -1.0, 1.0))  # Ensure the value is within the valid range for acos