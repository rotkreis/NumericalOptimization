import numpy as np
import scipy.optimize
from leastsquares import GN, LMF, Dogleg
u = np.array([
4.0e0, 2.0e0, 1.0e0, 5.0e-1, 2.5e-1, 1.67e-1, 1.25e-1,
1.0e-1, 8.33e-2, 7.14e-2, 6.25e-2
])
y = np.array([1.957e-1, 1.947e-1, 1.735e-1, 1.6e-1, 8.44e-2, 6.27e-2,
4.56e-2, 3.42e-2, 3.23e-2, 2.35e-2, 2.46e-2])
def r(x):
    res = np.zeros(11)
    for i in range(0,11):
        res[i] = y[i] - x[0]*(u[i]**2+u[i]*x[1])/(u[i]**2+u[i]*x[2]+x[3])
    return res
def jac(x):
    res = np.zeros((11,4))
    for i in range(0,11):
        denom = (u[i]**2+u[i]*x[2]+x[3])
        nom = x[0]*(u[i]**2+u[i]*x[1])
        res[i][0] = -(u[i]**2+u[i]*x[1]) / denom
        res[i][1] = -x[0]*u[i] / denom
        res[i][2] = u[i] * nom / denom**2
        res[i][3] = nom / denom ** 2
    return res

xs = np.array([2.5e-1, 3.9e-1, 4.15e-1, 3.9e-1])
xs = 10*xs
#GN(r, xs, jac, search = True)
#LMF(r, xs, jac)
Dogleg(r, xs, jac, maxiter = 2000)
Dogleg(r, xs, jac, maxiter = 4000)
Dogleg(r, xs, jac, maxiter = 6000)
res = scipy.optimize.leastsq(r, xs, Dfun = jac)
print res
x = res[0]
x1 = x
#print 0.5*np.dot(r(x1),r(x1))
