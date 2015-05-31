import numpy as np
import scipy.optimize
from leastsquares import GN, LMF, Dogleg
from largeresidual import DGW

def brown(n):
    def r(x):
        res = np.zeros(n)
        for i in range(0,n):
            i = i + 1
            [x1,x2,x3,x4] = x
            temp = i / 5.0
            res[i - 1] = (x1 + temp*x2 - np.exp(temp))**2 + (x3 + np.sin(temp)*x4 - np.cos(temp))**2
        return res
    return r

def brown_jac(n):
    def jac(x):
        res = np.zeros((n,4))
        for i in range(0,n):
            i = i + 1
            [x1,x2,x3,x4] = x
            temp= i / 5.0
            tmp1 = x1 + temp*x2 - np.exp(temp)
            tmp2 = x3 + np.sin(temp)*x4 - np.cos(temp)
            res[i-1][0] = 2 * tmp1
            res[i-1][1] = 2 * tmp1 * temp
            res[i-1][2] = 2 * tmp2
            res[i-1][3] = 2 * np.sin(temp) * tmp2
        return res
    return jac

xs = np.array([25,5,-5,1])

n = 10
r = brown(n)
jac = brown_jac(n)

DGW(r, xs, jac)
GN(r, xs, jac)
LMF(r, xs, jac)
res = scipy.optimize.leastsq(r, xs, Dfun = jac)
print res

