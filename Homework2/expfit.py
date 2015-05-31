import numpy as np
import scipy.optimize
from leastsquares import GN, LMF, Dogleg
from largeresidual import DGW
y = np.array([
8.44e-1, 9.08e-1, 9.32e-1, 9.36e-1, 9.25e-1, 9.08e-1,
8.81e-1, 8.5e-1, 8.18e-1, 7.84e-1, 7.51e-1, 7.18e-1, 6.85e-1,
6.58e-1, 6.28e-1, 6.03e-1, 5.8e-1, 5.58e-1, 5.38e-1, 5.22e-1,
5.06e-1, 4.9e-1, 4.78e-1, 4.67e-1, 4.57e-1, 4.48e-1, 4.38e-1,
4.31e-1, 4.24e-1, 4.2e-1, 4.14e-1, 4.11e-1, 4.06e-1])

def r(x):
    res = np.zeros(33)
    for i in range(0,33):
        ti = 10 * (i)
        res[i] = y[i] - (x[0] + x[1]*np.exp(-ti*x[3]) + x[2]*np.exp(-ti*x[4]))
        if np.isinf(res[i]):
            res[i] = 1e12
    return res
def jac(x):
    res = np.zeros((33,5))
    for i in range(0,33):
        ti = 10 * i
        res[i][0] = -1
        res[i][1] = -np.exp(-ti*x[3])
        res[i][2] = -np.exp(-ti*x[4])
        res[i][3] = x[1]*ti*np.exp(-ti*x[3])
        res[i][4] = x[2]*ti*np.exp(-ti*x[4])
    for i in range(0,33):
        for j in range(0,5):
            if np.isinf(res[i][j]):
                res[i][j] = 1e12
    return res
def h(x):
    res = np.zeros((5,5))
    [x1,x2,x3,x4,x5] = x
    for i in range(0,33):
        ti = 10 * i
        hess = np.zeros((5,5))
        hess[1][3] =  ti*np.exp(-ti*x4)
        hess[2][4] =  ti*np.exp(-ti*x5)
        hess[3][1] =  ti*np.exp(-ti*x4)
        hess[3][3] =  -ti**2*x2*np.exp(-ti*x4)
        hess[4][2] =  ti*np.exp(-ti*x5)
        hess[4][4] =  -ti**2*x3*np.exp(-ti*x5)
        res += hess * r(x)[i]
    return res
xs = np.array([5.0e-1, 1.5, -1, 1.0e-2, 2.0e-2])
GN(r, xs, jac, search = True, ave = 1e-6)
#LMF(r, xs, jac)
#DGW(r, xs, jac)
DGW(r, xs, jac)
DGW(r, xs, jac, h, method = "BFGS")
#Dogleg(r, xs, jac)
res = scipy.optimize.leastsq(r, xs, Dfun = jac)
print res
#print h(xs)
