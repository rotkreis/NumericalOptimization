import numpy as np
import scipy.optimize
from leastsquares import GN, LMF, Dogleg
from largeresidual import DGW

def sampson(n):
    def r(x):
        res = np.zeros(n)
        for i in range(0,n):
            j = i + 1
            res[i] = 2 + 2 * j - np.exp(j*x[0]) - np.exp(j*x[1])
        return res
    return r

def sampson_der(n):
    def jac(x):
        res = np.zeros((n,2))
        for i in range(0,n):
            j = i + 1
            res[i][0] = -j * np.exp(j * x[0])
            res[i][1] = -j * np.exp(j * x[1])
        return res
    return jac

sol10 = np.array([0.2578, 0.2578])
xs = np.array([0.3, 0.4])

n = 10
r = sampson(10)
jac = sampson_der(10)

#GN(r, xs, jac)
DGW(r, xs, jac)
