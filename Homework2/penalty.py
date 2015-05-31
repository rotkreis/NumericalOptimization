"""
    Methods for Constrained Problems
    Two kinds of Penalty Functions
    Specially for Least Squares
"""
import numpy as np
from scipy import linalg as LA
import scipy.optimize
from leastsquares import GN
def BoundPenalty(r,x0,jac,cons,u0=1,e0=1e-5,search=True,ave=1e-5, maxiter=1000, diag=False):
    """
    BoundPenalty solver
    This comment is rather meaningless
    """
    def f(x):
        return .5 * np.dot(r(x),r(x))
    def fprime(x):
        return np.dot(np.transpose(jac(x)),r(x))
    def penaly(x,u):
        sum = 0
        n = len(cons(x))
        for i in range(0, n):
            sum += 1.0 / cons(x)[i]
        return u * sum
    def b(x,u):
        return f(x) + penaly(x,u)
    xk = x0
    uk = u0
    while penaly(xk, uk) > e0 and count <= maxiter:
       xk = GN(dd)
