import numpy as np
import scipy.optimize
from leastsquares import GN, LMF, Dogleg
from largeresidual import DGW

def box(n):
    def r(x):
        res = np.zeros(n)
        for i in range(0,n):
            res[i] = ( np.exp(-(i+1)*x[0]/10.0) - np.exp(-(i+1) *x[1] /10.0) +
                      (np.exp(-(i+1))-np.exp(-(i+1)/10)) * x[2] )
        return res
    return r

def box_jac(n):
    def jac(x):
        res = np.zeros((n,3))
        for i in range(0,n):
            j = i + 1
            temp = j / 10.0
            res[i][0] = -temp * np.exp(-temp * x[0])
            res[i][1] = temp * np.exp(-temp * x[1])
            res[i][2] = np.exp(-j) - np.exp(-temp)
        return res
    return jac
sol = np.array([1.0, 10.0, 1.0])
xs = np.array([0, 10.0, 20.0])
n = 50
#xs = sol
r = box(n)
jac = box_jac(n)
def f(x):
    return .5 * np.dot(r(x),r(x))
def fprime(x):
    return np.dot(np.transpose(jac(x)),r(x))

#GN(r, xs, jac)
#LMF(r,xs,jac)
#Dogleg(r, xs, jac)
#DGW(r, xs, jac)
#DGW(r, xs, jac, method = "Bigg")
DGW(r, xs, jac, method = "BFGS", maxiter = 1000)
#res = scipy.optimize.leastsq(r, xs, Dfun = jac)
#print res[0]
#print f(res[0])
#print f(sol)
res = scipy.optimize.minimize(f, xs, method='BFGS',jac = fprime,
               options={'disp':True})
print res.x

