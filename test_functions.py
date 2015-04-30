#from scipy.optimize.linesearch import line_search_wolfe2 as line_search
#import numpy as np
#from scipy import linalg as LA
#import scipy.optimize
#from watson import watson, watson_der, watson_hess
from sympy import diff, Symbol, symbols
#import sympy as sp

#x1,x2,x3,x4,x5 = symbols('x1 x2 x3 x4 x5')
#r1,r2,r3,r4,r4,r5,r6,r7,r8,r9,r10 = symbols('r1 r2 r3 r4 r4 r5 r6 r7 r8 r9 r10')
#r = Symbol('r')
#f1 = x1*x2 + x1 - 3*x5
#f2 = 2*x1*x2 + x1 + 2*r10*x2**2 + x2*x3**2 + r7*x2*x3 + r9*x2*x4 + r8*x2 - r*x5
#f3 = 2*x2*x3**2 + r7*x2*x3 + 2*r5*x3**2 + r6*x3 - 8*x5
#f4 = r9*x2*x4 + 2*x4**2 - 4*r*x5
#f5 = (x1*x2 + x1 + r10*x2**2 + x2*x3**2 + r7*x2*x3 + r9*x2*x4 + r8*x2 +
      #r5*x3**2 + r6*x3 + x4**2 - 1)

#f = f1**2 + f2**2 + f3**2 + f4**2 + f5**2
##for x in [x1,x2,x3,x4,x5]:
    ##print diff(f,x)

#result = open('res_propane','w')
#x_array = [x1,x2,x3,x4,x5]
##for x in x_array:
    ##print 'der[' + str(x_array.index(x)) + '] = ',
    ##print diff(f,x);
#for x in x_array:
    #for y in x_array:
        #print ('hess['+str(x_array.index(x)) +'][' +
               #str(x_array.index(y)) + '] = '),
        #print diff(f,x,y)
#t = diff(f,x1,x2)
## in vim use :%s/koko\n//g
#print 'a[1] = koko'
#print t

#dipole 2.6
x1,x2,x3,x4,x5,x6,x7,x8 = symbols('x1 x2 x3 x4 x5 x6 x7 x8')
sig_x,sig_y, sig_a,sig_b,sig_c,sig_d,sig_e,sig_f = symbols('sig_x sig_y sig_a sig_b sig_c sig_d sig_e sig_f')
f1 = x1 + x2 - sig_x
f2 = x3 + x4 - sig_y
f3 = x5*x1 + x6*x2 - x7*x3 - x8*x4 - sig_a
f4 = x7*x1 + x8*x2 + x5*x3 + x6*x4 - sig_b
f5 = x1*(x5**2-x7**2) - 2*x3*x5*x7 + x2*(x6**2-x8**2) - 2*x4*x6*x8 - sig_c
f6 = x3*(x5**2-x7**2) + 2*x1*x5*x7 + x4*(x6**2-x8**2) + 2*x2*x6*x8 - sig_d
f7 = x1*x5*(x5**2 - 3*x7**2) + x3*x7*(x7**2 - 3*x5**2) + x2*x6*(x6**2 - 3*x8**2) + x4*x8*(x8**2 - 3*x6**2) - sig_e
f8 = x3*x5*(x5**2 - 3*x7**2) - x1*x7*(x7**2 - 3*x5**2) + x4*x6*(x6**2 - 3*x8**2) - x2*x8*(x8**2 - 3*x6**2) - sig_f

f = f1**2 + f2**2 + f3**2 + f4**2 + f5**2 + f6**2 + f7**2 + f8**2
result = open('res_dipole','w')
x_array = [x1,x2,x3,x4,x5,x6,x7,x8]
for i,x in enumerate(x_array):
    print 'der[' + str(i) + '] = ',
    print diff(f,x)
for i,x in enumerate(x_array):
    for j,y in enumerate(x_array):
        print 'hess[' + str(i) + '][' + str(j) + '] = ',
        print diff(f,x,y)




















