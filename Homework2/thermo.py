import numpy as np
import scipy.optimize
from leastsquares import GN, LMF, Dogleg

y = np.array([3.478e4, 2.861e4, 2.365e4, 1.963e4, 1.637e4, 1.372e4,
1.154e4, 9.744e3, 8.261e3, 7.03e3, 6.005e3, 5.147e3, 4.427e3,
3.82e3, 3.307e3, 2.872e3])

def r(x):
    res = np.zeros(16)
    for i in range(0,16):
        ti = 5 + 45*(i+1)
        res[i] = x[0]*np.exp(x[1] / (ti+x[2])) - y[i]
    return res
def jac(x):
    res = np.zeros((16,3))
    for i in range(0,16):
        ti = 5 + 45*(i+1)
        e = np.exp(x[1] / (ti + x[2]))
        res[i][0] = e
        res[i][1] = x[0] / (ti+x[2]) * e
        res[i][2] = -x[0]*x[1] / (ti+x[2])**2 * e
    return res

xs = np.array([2.0e-2, 4.0e3, 2.5e2])
x0 = xs

#GN(r, xs, jac, search = True)
#LMF(r, xs, jac)
Dogleg(r, xs, jac)
res = scipy.optimize.leastsq(r, xs, Dfun = jac)
print res
