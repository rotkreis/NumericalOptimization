import numpy as np
import scipy.optimize
from cgmethods import cgsolve
from cgmethods import BBsolve
def s303(n):
    def f(x):
        sum = 0
        s1= 0
        for i in range(0, len(x)):
            sum += x[i]**2
            s1 += (i+1)*x[i]/2.0
        sum += s1**2
        sum += s1**4
        return sum
    return f
def s303_der(n):
    def f(x):
        der = np.zeros(n)
        sum = 0
        for i in range(0, n):
            sum += (i+1)*x[i]/2.0
        for i in range(0, n):
            der[i] = 2*x[i] + sum*(i+1) + 2*sum**3*(i+1)
        return der
    return f
def gencube(n):
    def f(x):
        sum = 0
        sum += (x[0]-1)**2
        for i in range(1,len(x)):
            sum += 100 * (x[i]-x[i-1]**3)**2
        return sum
    return f
def gencube_der(n):
    def f(x):
        der = np.zeros(n)
        der[0] = 2*(x[0]-1) - 200 * (x[1]-x[0]**3) * 3*x[0]**2
        for i in range(1,n - 1):
            der[i] = 200*((x[i] - x[i-1]**3) - (x[i+1] - x[i]**3)*3*x[i]**2)
        der[n-1] = 200*(x[n-1] - x[n-2]**3)
        return der
    return f
def gensin(n):
    def f(x):
        sum = 0
        c1 = 1e-4
        c2 = 4
        for i in range(0, n-1):
            sum += (x[i+1]-np.sin(x[i]))**2 / c1 + x[i]**2 / c2
        return sum
    return f
def gensin_der(n):
    def f(x):
        c1 = 1e-4
        c2 = 4
        der = np.zeros(n)
        der[0] = 2*(x[1]-np.sin(x[0]))/c1*(-np.cos(x[0])) + 2*x[0]/c2
        for i in range(1,n-1):
            der[i] = (2*((x[i+1]-np.sin(x[i]))/c1*(-np.cos(x[i])) + x[i]/c2 +
                      (x[i]-np.sin(x[i-1]))/c1))
        der[n-1] = 2*(x[n-1] - np.sin(x[n-2]))/c1
        return der
    return f

n = 81
#x0 =0.8 * np.ones(n)
x0 = np.ones(n) * 0.1
#print s303_der(n)(x0)
#cgsolve(gensin(n),x0,gensin_der(n), "FR")
#cgsolve(gensin(n),x0,gensin_der(n), "PRP")
#cgsolve(gensin(n),x0,gensin_der(n), "PRP+")
#cgsolve(gensin(n),x0,gensin_der(n), "CD")
#cgsolve(gensin(n),x0,gensin_der(n), "")
#BBsolve(gencube(n),x0,gencube_der(n), method = 2)
#BBsolve(s303(n),x0, s303_der(n), method = 1)
cgsolve(s303(n),x0, s303_der(n), method = "FR")
#cgsolve(gencube(n),x0,gencube_der(n), "FR")
#cgsolve(gencube(n),x0,gencube_der(n), "PRP")
#cgsolve(gencube(n),x0,gencube_der(n), "PRP+")
#cgsolve(gencube(n),x0,gencube_der(n), "CD")
#cgsolve(gencube(n),x0,gencube_der(n), "DY")
res = scipy.optimize.minimize(s303(n), x0, method='BFGS',jac = s303_der(n),
               options={'disp':True})
print res.x

#print scipy.optimize.check_grad(gencube(n), gencube_der(n), x0)
#print scipy.optimize.approx_fprime(x0,gencube(n),1e-5)
#der_n = gencube_der(n)
#print der_n(x0)
