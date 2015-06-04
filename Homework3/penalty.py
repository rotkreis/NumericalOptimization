"""
    Methods for Constrained Problems
    Two kinds of Penalty Functions
    Specially for Least Squares
"""
import numpy as np
from scipy import linalg as LA
import scipy.optimize
from NewtonMethods import BFGS
def print_res(msg, iter, totalfc, totalgc,  gk, xk, fk):
    print msg
    print "Iterations: ",
    print iter
    print "Function Evaluations: ",
    print totalfc
    print "Gradient Evaluations: ",
    print totalgc
    print 'gk = '
    print gk
    print 'xk = '
    print xk
    print 'fk = '
    #print '%.6f' % fk
    print fk
    print " "


def BoundPenalty(f,x0,fprime,cons, cons_der,
                 u0=1,e0=1e-5,search=True,ave=1e-8,
                 maxiter=1000, diag=False, disp = True):
    """
    BoundPenalty solver
    This comment is rather meaningless
    """
    def penaly(x,u):
        sum = 0
        n = len(cons(x))
        for i in range(0, n):
            sum += 1.0 / cons(x)[i]
        return u * sum
    def b(u):
        def res(x):
            return f(x) + penaly(x,u)
        return res
    def bprime(u):
        def res(x):
            return fprime(x) + u * cons_der(x)
        return res
    xk = x0
    uk = u0
    count = 0
    while penaly(xk, uk) > ave and count <= maxiter:
        res = BFGS(b(uk), xk, bprime(uk),  ave = ave, disp = False)
        xk = res[0]
        warnflag = res[1]
        msg = res[2]
        if warnflag != 0:
            print msg
            #break
        uk = uk / 10.0
        count += 1
        print uk

    return xk,uk, penaly(xk,uk), count, f(xk), fprime(xk)

