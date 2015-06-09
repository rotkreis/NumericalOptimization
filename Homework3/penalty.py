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
    def penalty(x,u):
        sum = 0
        n = len(cons(x))
        for i in range(0, n):
            temp = cons(x)[i]
            if temp != 0:
                sum += 1.0 / temp
        return u * sum
    def b(u):
        def res(x):
            return f(x) + penalty(x,u)
        return res
    def bprime(u):
        def res(x):
            return fprime(x) + u * cons_der(x)
        return res
    xk = x0
    uk = u0
    count = 0
    while penalty(xk, uk) > ave and count <= maxiter:
        #return iter,totalfc, totalgc, fpk,xk, f(xk), warnflag, msg
        res = BFGS(b(uk), xk, bprime(uk),  ave = ave, disp = False)
        xk = res[-4]
        warnflag = res[-2]
        msg = res[-1]
        if warnflag != 0:
            print msg
            #break
        uk = uk / 10.0
        count += 1
        print uk
    if disp == True:
        print_res("BoundPenalty fin", res[0],res[1],res[2],res[3],res[4],res[5])
    return xk,uk, penalty(xk,uk), count, f(xk), fprime(xk)


def Lagrange(f,x0,fprime,ins,ins_der,sigma0,gamma0,
             e0=1e-5,ave=1e-8,maxiter=1000,diag=False,disp=True):
    "Implement only inequalities"
    def l(gamma,sigma):
        def temp(x):
            return f(x) + penalty(x,gamma,sigma)
        return temp
    def penalty(x,gamma,sigma):
        return None # TEMP!
    #Return FORMAT: iter,totalfc, totalgc, fpk,xk, f(xk), warnflag, msg
    while LA.norm(ins(xk)) > ave and count <= maxiter:
        res = BFGS(l(gamma_k,sigma_k), xk, lprime(gamma_k,sigma_k),ave=e0,disp=False)
        xk = res[-4]
        warnflag = res[-2]
        msg = res[-1]
        if warnflag != 0:
            print msg
        #update gamma_k
        sigma_k = 10 * sigma_k
        count += 1
    if disp == True:
        print_res("Lagrange fin", res[0],res[1],res[2],res[3],res[4],res[5])
    return xk, sigma_k,gamma_k,  LA.norm(ins(xk)), count, f(xk), fprime(xk)


def ExPenalty(f,x0,fprime,eqns,eqns_der,ins,ins_der,
              sigma0=1,e0=1e-5,ave=1e-8,maxiter=1000,
              diag=False, disp=True):
    def penalty(x,sigma):
        sum = 0
        n = len(eqns(x))
        for i in range(0,n):
            sum += (eqns(x)[i])**2
            sum += np.min([ins(x)[i],0])**2
        return .5 * sigma * sum
    def p(sigma):
        def temp(x):
            return f(x) + penalty(x,sigma)
        return temp
    def pprime(sigma):
        def temp(x):
            return fprime(x) + 0.5 * sigma * (eqns_der(x)+ins_der(x))
        return temp
    xk = x0
    sigma_k = sigma0
    count = 0
    while LA.norm(eqns(xk)) > ave and count <= maxiter:
        #return iter,totalfc, totalgc, fpk,xk, f(xk), warnflag, msg
        res = BFGS(p(sigma_k), xk, pprime(sigma_k),ave=ave,disp=False)
        xk = res[-4]
        warnflag = res[-2]
        msg = res[-1]
        if warnflag != 0:
            print msg
            #break
        sigma_k = 10 * sigma_k
        count += 1
        print sigma_k
    if disp == True:
        print_res("Exterior fin", res[0],res[1],res[2],res[3],res[4],res[5])
    return xk, sigma_k, LA.norm(eqns(xk)), count, f(xk), fprime(xk)
