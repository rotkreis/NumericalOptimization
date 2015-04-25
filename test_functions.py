#from scipy.optimize.linesearch import line_search_wolfe2 as line_search
#import numpy as np
#from scipy import linalg as LA
#import scipy.optimize
#from watson import watson, watson_der, watson_hess
from sympy import *

x1,x2,x3,x4,x5 = symbols('x1 x2 x3 x4 x5')
r1,r2,r3,r4,r4,r5,r6,r7,r8,r9,r10 = symbols('r1 r2 r3 r4 r4 r5 r6 r7 r8 r9 r10')
r = Symbol('r')
f1 = x1*x2 + x1 - 3*x5
f2 = 2*x1*x2 + x1 + 2*r10*x2**2 + x2*x3**2 + r7*x2*x3 + r9*x2*x4 + r8*x2 - r*x5
f3 = 2*x2*x3**2 + r7*x2*x3 + 2*r5*x3**2 + r6*x3 - 8*x5
f4 = r9*x2*x4 + 2*x4**2 - 4*r*x5
f5 = (x1*x2 + x1 + r10*x2**2 + x2*x3**2 + r7*x2*x3 + r9*x2*x4 + r8*r2 +
      r5*r3**2 + r6*r3 + x4**2 - 1)

f = f1**2 + f2**2 + f3**2 + f4**2 + f5**2
for x in [x1,x2,x3,x4,x5]:
    print diff(f,x)
for x in [x1,x2,x3,x4,x5]:
    for y in [x1,x2,x3,x4,x5]:
        print diff(f,x,y);
def propane(x1,x2,x3,x4,x5):
    return f
print propane(1,2,3,4,5);
