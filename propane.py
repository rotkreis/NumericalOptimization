import numpy as np
k5 = 1.930e-1
k6 = 2.597e-3
k7 = 3.448e-3
k8 = 1.799e-5
k9 = 2.155e-4
k10 = 3.846e-5
p = 40
r = 10
r5 = k5
r6 = k6 * p ** (-0.5)
r7 = k7 * p ** (-0.5)
r8 = k8 / p
r9 = k9 * p ** (-0.5)
r10 = k10 / p
def propane(x):
    x1 = x[0]
    x2 = x[1]
    x3 = x[2]
    x4 = x[3]
    x5 = x[4]
    f1 = x1*x2 + x1 - 3*x5
    f2 = 2*x1*x2 + x1 + 2*r10*x2**2 + x2*x3**2 + r7*x2*x3 + r9*x2*x4 + r8*x2 - r*x5
    f3 = 2*x2*x3**2 + r7*x2*x3 + 2*r5*x3**2 + r6*x3 - 8*x5
    f4 = r9*x2*x4 + 2*x4**2 - 4*r*x5
    f5 = (x1*x2 + x1 + r10*x2**2 + x2*x3**2 + r7*x2*x3 + r9*x2*x4 + r8*x2 +
          r5*x3**2 + r6*x3 + x4**2 - 1)
    return f1**2 + f2**2 + f3**2 + f4**2 + f5**2

def propane_der(x):
    x1 = x[0]
    x2 = x[1]
    x3 = x[2]
    x4 = x[3]
    x5 = x[4]
    der = np.zeros(5)
    der[0] =  (2*x2 + 2)*(x1*x2 + x1 - 3*x5) + (2*x2 + 2)*(r10*x2**2 + r5*x3**2 + r6*x3 + r7*x2*x3 + r8*x2 + r9*x2*x4 + x1*x2 + x1 + x2*x3**2 + x4**2 - 1) + (4*x2 + 2)*(-r*x5 + 2*r10*x2**2 + r7*x2*x3 + r8*x2 + r9*x2*x4 + 2*x1*x2 + x1 + x2*x3**2)
    der[1] =  2*r9*x4*(-4*r*x5 + r9*x2*x4 + 2*x4**2) + 2*x1*(x1*x2 + x1 - 3*x5) + (2*r7*x3 + 4*x3**2)*(2*r5*x3**2 + r6*x3 + r7*x2*x3 + 2*x2*x3**2 - 8*x5) + (4*r10*x2 + 2*r7*x3 + 2*r8 + 2*r9*x4 + 2*x1 + 2*x3**2)*(r10*x2**2 + r5*x3**2 + r6*x3 + r7*x2*x3 + r8*x2 + r9*x2*x4 + x1*x2 + x1 + x2*x3**2 + x4**2 - 1) + (8*r10*x2 + 2*r7*x3 + 2*r8 + 2*r9*x4 + 4*x1 + 2*x3**2)*(-r*x5 + 2*r10*x2**2 + r7*x2*x3 + r8*x2 + r9*x2*x4 + 2*x1*x2 + x1 + x2*x3**2)
    der[2] =  (2*r7*x2 + 4*x2*x3)*(-r*x5 + 2*r10*x2**2 + r7*x2*x3 + r8*x2 + r9*x2*x4 + 2*x1*x2 + x1 + x2*x3**2) + (4*r5*x3 + 2*r6 + 2*r7*x2 + 4*x2*x3)*(r10*x2**2 + r5*x3**2 + r6*x3 + r7*x2*x3 + r8*x2 + r9*x2*x4 + x1*x2 + x1 + x2*x3**2 + x4**2 - 1) + (8*r5*x3 + 2*r6 + 2*r7*x2 + 8*x2*x3)*(2*r5*x3**2 + r6*x3 + r7*x2*x3 + 2*x2*x3**2 - 8*x5)
    der[3] =  2*r9*x2*(-r*x5 + 2*r10*x2**2 + r7*x2*x3 + r8*x2 + r9*x2*x4 + 2*x1*x2 + x1 + x2*x3**2) + (2*r9*x2 + 4*x4)*(r10*x2**2 + r5*x3**2 + r6*x3 + r7*x2*x3 + r8*x2 + r9*x2*x4 + x1*x2 + x1 + x2*x3**2 + x4**2 - 1) + (2*r9*x2 + 8*x4)*(-4*r*x5 + r9*x2*x4 + 2*x4**2)
    der[4] =  -8*r*(-4*r*x5 + r9*x2*x4 + 2*x4**2) - 2*r*(-r*x5 + 2*r10*x2**2 + r7*x2*x3 + r8*x2 + r9*x2*x4 + 2*x1*x2 + x1 + x2*x3**2) - 32*r5*x3**2 - 16*r6*x3 - 16*r7*x2*x3 - 6*x1*x2 - 6*x1 - 32*x2*x3**2 + 146*x5
    return der

def propane_hess(x):
    x1 = x[0]
    x2 = x[1]
    x3 = x[2]
    x4 = x[3]
    x5 = x[4]
    hess = np.zeros((5,5))
    hess[0][0] =  2*(2*(x2 + 1)**2 + (2*x2 + 1)**2)
    hess[0][1] =  2*(-2*r*x5 + 5*r10*x2**2 + r5*x3**2 + r6*x3 + 3*r7*x2*x3 + 3*r8*x2 + 3*r9*x2*x4 + 6*x1*x2 + x1*(x2 + 1) + 4*x1 + 3*x2*x3**2 + x4**2 - 3*x5 + (x2 + 1)*(2*r10*x2 + r7*x3 + r8 + r9*x4 + x1 + x3**2) + (2*x2 + 1)*(4*r10*x2 + r7*x3 + r8 + r9*x4 + 2*x1 + x3**2) - 1)
    hess[0][2] =  2*(x2*(r7 + 2*x3)*(2*x2 + 1) + (x2 + 1)*(2*r5*x3 + r6 + r7*x2 + 2*x2*x3))
    hess[0][3] =  2*(r9*x2*(2*x2 + 1) + (x2 + 1)*(r9*x2 + 2*x4))
    hess[0][4] =  -2*(r*(2*x2 + 1) + 3*x2 + 3)
    hess[1][0] =  2*(-2*r*x5 + 5*r10*x2**2 + r5*x3**2 + r6*x3 + 3*r7*x2*x3 + 3*r8*x2 + 3*r9*x2*x4 + 6*x1*x2 + x1*(x2 + 1) + 4*x1 + 3*x2*x3**2 + x4**2 - 3*x5 + (x2 + 1)*(2*r10*x2 + r7*x3 + r8 + r9*x4 + x1 + x3**2) + (2*x2 + 1)*(4*r10*x2 + r7*x3 + r8 + r9*x4 + 2*x1 + x3**2) - 1)
    hess[1][1] =  2*(4*r10*(-r*x5 + 2*r10*x2**2 + r7*x2*x3 + r8*x2 + r9*x2*x4 + 2*x1*x2 + x1 + x2*x3**2) + 2*r10*(r10*x2**2 + r5*x3**2 + r6*x3 + r7*x2*x3 + r8*x2 + r9*x2*x4 + x1*x2 + x1 + x2*x3**2 + x4**2 - 1) + r9**2*x4**2 + x1**2 + x3**2*(r7 + 2*x3)**2 + (2*r10*x2 + r7*x3 + r8 + r9*x4 + x1 + x3**2)**2 + (4*r10*x2 + r7*x3 + r8 + r9*x4 + 2*x1 + x3**2)**2)
    hess[1][2] =  2*(x2*(r7 + 2*x3)*(4*r10*x2 + r7*x3 + r8 + r9*x4 + 2*x1 + x3**2) + x3*(r7 + 2*x3)*(4*r5*x3 + r6 + r7*x2 + 4*x2*x3) + (r7 + 2*x3)*(-r*x5 + 2*r10*x2**2 + r7*x2*x3 + r8*x2 + r9*x2*x4 + 2*x1*x2 + x1 + x2*x3**2) + (r7 + 2*x3)*(r10*x2**2 + r5*x3**2 + r6*x3 + r7*x2*x3 + r8*x2 + r9*x2*x4 + x1*x2 + x1 + x2*x3**2 + x4**2 - 1) + (r7 + 4*x3)*(2*r5*x3**2 + r6*x3 + r7*x2*x3 + 2*x2*x3**2 - 8*x5) + (2*r5*x3 + r6 + r7*x2 + 2*x2*x3)*(2*r10*x2 + r7*x3 + r8 + r9*x4 + x1 + x3**2))
    hess[1][3] =  2*(r9*x2*(4*r10*x2 + r7*x3 + r8 + r9*x4 + 2*x1 + x3**2) + r9*x4*(r9*x2 + 4*x4) + r9*(-4*r*x5 + r9*x2*x4 + 2*x4**2) + r9*(-r*x5 + 2*r10*x2**2 + r7*x2*x3 + r8*x2 + r9*x2*x4 + 2*x1*x2 + x1 + x2*x3**2) + r9*(r10*x2**2 + r5*x3**2 + r6*x3 + r7*x2*x3 + r8*x2 + r9*x2*x4 + x1*x2 + x1 + x2*x3**2 + x4**2 - 1) + (r9*x2 + 2*x4)*(2*r10*x2 + r7*x3 + r8 + r9*x4 + x1 + x3**2))
    hess[1][4] =  -2*(4*r*r9*x4 + r*(4*r10*x2 + r7*x3 + r8 + r9*x4 + 2*x1 + x3**2) + 8*r7*x3 + 3*x1 + 16*x3**2)
    hess[2][0] =  2*(x2*(r7 + 2*x3)*(2*x2 + 1) + (x2 + 1)*(2*r5*x3 + r6 + r7*x2 + 2*x2*x3))
    hess[2][1] =  2*(x2*(r7 + 2*x3)*(4*r10*x2 + r7*x3 + r8 + r9*x4 + 2*x1 + x3**2) + x3*(r7 + 2*x3)*(4*r5*x3 + r6 + r7*x2 + 4*x2*x3) + (r7 + 2*x3)*(-r*x5 + 2*r10*x2**2 + r7*x2*x3 + r8*x2 + r9*x2*x4 + 2*x1*x2 + x1 + x2*x3**2) + (r7 + 2*x3)*(r10*x2**2 + r5*x3**2 + r6*x3 + r7*x2*x3 + r8*x2 + r9*x2*x4 + x1*x2 + x1 + x2*x3**2 + x4**2 - 1) + (r7 + 4*x3)*(2*r5*x3**2 + r6*x3 + r7*x2*x3 + 2*x2*x3**2 - 8*x5) + (2*r5*x3 + r6 + r7*x2 + 2*x2*x3)*(2*r10*x2 + r7*x3 + r8 + r9*x4 + x1 + x3**2))
    hess[2][2] =  2*(x2**2*(r7 + 2*x3)**2 + 2*x2*(-r*x5 + 2*r10*x2**2 + r7*x2*x3 + r8*x2 + r9*x2*x4 + 2*x1*x2 + x1 + x2*x3**2) + 4*(r5 + x2)*(2*r5*x3**2 + r6*x3 + r7*x2*x3 + 2*x2*x3**2 - 8*x5) + 2*(r5 + x2)*(r10*x2**2 + r5*x3**2 + r6*x3 + r7*x2*x3 + r8*x2 + r9*x2*x4 + x1*x2 + x1 + x2*x3**2 + x4**2 - 1) + (2*r5*x3 + r6 + r7*x2 + 2*x2*x3)**2 + (4*r5*x3 + r6 + r7*x2 + 4*x2*x3)**2)
    hess[2][3] =  2*(r9*x2**2*(r7 + 2*x3) + (r9*x2 + 2*x4)*(2*r5*x3 + r6 + r7*x2 + 2*x2*x3))
    hess[2][4] =  -2*(r*x2*(r7 + 2*x3) + 32*r5*x3 + 8*r6 + 8*r7*x2 + 32*x2*x3)
    hess[3][0] =  2*(r9*x2*(2*x2 + 1) + (x2 + 1)*(r9*x2 + 2*x4))
    hess[3][1] =  2*(r9*x2*(4*r10*x2 + r7*x3 + r8 + r9*x4 + 2*x1 + x3**2) + r9*x4*(r9*x2 + 4*x4) + r9*(-4*r*x5 + r9*x2*x4 + 2*x4**2) + r9*(-r*x5 + 2*r10*x2**2 + r7*x2*x3 + r8*x2 + r9*x2*x4 + 2*x1*x2 + x1 + x2*x3**2) + r9*(r10*x2**2 + r5*x3**2 + r6*x3 + r7*x2*x3 + r8*x2 + r9*x2*x4 + x1*x2 + x1 + x2*x3**2 + x4**2 - 1) + (r9*x2 + 2*x4)*(2*r10*x2 + r7*x3 + r8 + r9*x4 + x1 + x3**2))
    hess[3][2] =  2*(r9*x2**2*(r7 + 2*x3) + (r9*x2 + 2*x4)*(2*r5*x3 + r6 + r7*x2 + 2*x2*x3))
    hess[3][3] =  2*(-16*r*x5 + 2*r10*x2**2 + 2*r5*x3**2 + 2*r6*x3 + 2*r7*x2*x3 + 2*r8*x2 + r9**2*x2**2 + 6*r9*x2*x4 + 2*x1*x2 + 2*x1 + 2*x2*x3**2 + 10*x4**2 + (r9*x2 + 2*x4)**2 + (r9*x2 + 4*x4)**2 - 2)
    hess[3][4] =  -2*r*(5*r9*x2 + 16*x4)
    hess[4][0] =  -2*(r*(2*x2 + 1) + 3*x2 + 3)
    hess[4][1] =  -2*(4*r*r9*x4 + r*(4*r10*x2 + r7*x3 + r8 + r9*x4 + 2*x1 + x3**2) + 8*r7*x3 + 3*x1 + 16*x3**2)
    hess[4][2] =  -2*(r*x2*(r7 + 2*x3) + 32*r5*x3 + 8*r6 + 8*r7*x2 + 32*x2*x3)
    hess[4][3] =  -2*r*(5*r9*x2 + 16*x4)
    hess[4][4] =  2*(17*r**2 + 73)
    return hess
