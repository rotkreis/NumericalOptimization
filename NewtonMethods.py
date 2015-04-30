from scipy.optimize.linesearch import line_search_BFGS,line_search_armijo, line_search_wolfe1
from scipy.optimize.linesearch import line_search_wolfe2
from scipy.optimize import rosen, rosen_der, rosen_hess
import numpy as np
from scipy import linalg as LA
import scipy.optimize
from watson import watson, watson_der, watson_hess
from propane import propane, propane_der, propane_hess
from cluster import cluster, cluster_der
class _LineSearchError(RuntimeError):
    pass
def line_search_wolfe12(f, fprime, xk, pk, gfk, old_fval, old_old_fval,
                         **kwargs):
    """
    same as line_search_wolfe1, but fall back to line_search_wolfe2 if
    suitable step length is not found, and raise an exception if a
    suitable step length is not found.

    raises
    ------
    _linesearcherror
        if no suitable step size is found

    """
    ret = line_search_wolfe1(f, fprime, xk, pk, gfk,
                             old_fval, old_old_fval,
                             **kwargs)

    if ret[0] is None:
        # line search failed: try different one.
        ret = line_search_wolfe2(f, fprime, xk, pk, gfk,
                                 old_fval, old_old_fval)

    if ret[0] is None:
        raise _LineSearchError()

    return ret

def stepNewton(f, x0, fprime, fhess, ave = 1e-5, maxiter = 2000):
    """
    my simple minimization algorithm, testing
    using newton Method
    """
    iter = 0
    fpk = fprime(x0)
    xk = x0
    old_fval = f(x0)
    old_old_fval = None
    warnflag = 0
    while (LA.norm(fpk) > ave and iter < maxiter):
        hess = fhess(xk)
        dk = LA.solve(hess, fpk * -1)
        try:
            step, fc, gc, old_fval, old_old_fval, gfkp1 = \
                line_search_wolfe12(f, fprime, xk, dk, fpk, old_fval, old_old_fval)
        except _LineSearchError:
            warnflag = 2
            break
        xk += dk * step
        fpk = fprime(xk)
        iter += 1
    if warnflag == 2:
        print "Line search error in stepNewton at iteration: ",
        print iter
    else:
        print "StepNewton Finished"
        print "iterations:"
        print iter
        print 'gk = '
        print fpk
        print 'xk = '
        print xk

#stepNewton(watson, np.zeros(9), watson_der, watson_hess, 1e-8)

def adjustedNewton(f, x0, fprime, fhess, u = 1e-4, ave = 1e-6, maxiter = 2000):
    """
    my simple minimization algorithm, testing
    using LM - newton Method
    """
    iter = 0
    fpk = fprime(x0)
    xk = x0
    old_fval = f(x0)
    old_old_fval = None
    warnflag = 0
    while (LA.norm(fpk) > ave and iter < maxiter):
        hess = fhess(xk)
        while True:
            try:
                LA.cholesky(hess)
                break
            except LA.LinAlgError:
                hess += u * np.eye(len(x0),len(x0))
                u *= 2
        #while any(np.diag(l) < 0.01 ):
            #hess += u * np.eye(len(x0),len(x0))
            #u *= 2
        dk = LA.solve(hess, fpk * -1)
        try:
            step, fc, gc, old_fval, old_old_fval, gfkp1 = \
                line_search_wolfe12(f, fprime, xk, dk, fpk, old_fval, old_old_fval)
        except _LineSearchError:
            warnflag = 2
            break
        xk += dk * step
        iter += 1
        fpk = fprime(xk)
    if warnflag == 2:
        print "Line search error in adjustedNewton at iteration: ",
        print iter
    else:
        print "Finished"
        print "iterations:"
        print iter
        print 'gk = '
        print fpk
        print 'xk = '
        print xk
#adjustedNewton(watson, np.zeros(9), watson_der, watson_hess, 1e-8)

def SR1(f, x0, fprime, fhess, ave = 1e-6, maxiter = 1000):
    """
    Quasi-Newton: SR1
    fhess: inital hessian, converted to H0 for further compuatitions
    possible error: the update of sk, yk, hk
    """
    iter = 0
    fpk = fprime(x0)
    xk = x0
    old_fval = f(x0)
    old_old_fval = None
    I = np.eye(len(x0), dtype=int)
    hk = I
    warnflag = 0
    while (LA.norm(fpk) > ave and iter < maxiter):
        dk = - np.dot(hk,fpk)
        try:
            step, fc, gc, old_fval, old_old_fval, gfkp1 = \
                line_search_wolfe12(f, fprime, xk, dk, fpk, old_fval, old_old_fval)
        except _LineSearchError:
            warnflag = 2
            break
        sk = dk * step
        xk += sk
        yk = fprime(xk) - fpk
        fpk += yk
        temp = sk - np.dot(hk,yk)
        try:
            rhok = 1.0 / np.dot(temp, yk)
        except ZeroDivisionError:
            rhok = 1000.0
            print "Divided by Zero, SR1"
        if np.isinf(rhok):
            rhok = 1000.0
            print "Divided by Zero, SR1"
        hk += (np.outer(temp, temp) * rhok)
        iter += 1
    if warnflag == 2:
        print "Line search error in SR1 at iteration: ",
        print iter
    else:
        print "SR1 Finished"
        print "iterations:"
        print iter
        print 'gk = '
        print fpk
        print 'xk = '
        print xk


#print SR1(watson, np.ones(9), watson_der, watson_hess)

def BFGS(f, x0, fprime, fhess, ave = 1e-6, maxiter = 1000):
    """
    Quasi-Newton: BFGS
    fhess: inital hessian, converted to H0 for further compuatitions
    possible error: the update of sk, yk, hk
    """
    iter = 0
    fpk = fprime(x0)
    xk = x0
    old_fval = f(x0)
    #old_fval = None
    old_old_fval = None
    I = np.eye(len(x0), dtype=int)
    # Holy Shit! Initially inv(hessian!) now ok ! hahahhahaha
    hk = I
    warnflag = 0
    while (LA.norm(fpk) > ave and iter < maxiter):
        dk = - np.dot(hk,fpk)
        try:
            step, fc, gc, old_fval, old_old_fval, gfkp1 = \
                line_search_wolfe12(f, fprime, xk, dk, fpk, old_fval, old_old_fval)
        except _LineSearchError:
            warnflag = 2
            break
        xk1 = xk + step * dk
        sk = xk1 - xk
        xk = xk1
        fpk1 = fprime(xk1)
        yk = fpk1 - fpk
        fpk = fpk1
        try:
            rhok = 1.0 / np.dot(sk,yk)
        except ZeroDivisionError:
            rhok = 1000.0
            print "Divided by zero occured"
        if np.isinf(rhok):
            rhok = 1000.0
            print "Divided by zero occured"
        A1 = I - np.outer(sk, yk) * rhok
        A2 = I - np.outer(yk, sk) * rhok
        hk = np.dot(A1, np.dot(hk,A2)) + rhok * np.outer(sk,sk)
        iter += 1
    if warnflag == 2:
        print "Line search error in BFGS at iteration: ",
        print iter
    else:
        print "BFGS Finished"
        print "iterations:"
        print iter
        print 'gk = '
        print fpk
        print 'xk = '
        print xk

#BFGS(watson, np.zeros(9), watson_der, watson_hess)

def DFP(f, x0, fprime, fhess, ave = 1e-6, maxiter = 1000):
    """
    Quasi-Newton: DFP
    fhess: inital hessian, converted to H0 for further compuatitions
    possible error: the update of sk, yk, hk
    """
    iter = 0
    fpk = fprime(x0)
    print fpk
    xk = x0
    old_fval = f(x0)
    old_old_fval = None
    I = np.eye(len(x0), dtype=int)
    hk = I
    warnflag = 0
    while (LA.norm(fpk) > ave and iter < maxiter):
        dk = - np.dot(hk,fpk)
        try:
            step, fc, gc, old_fval, old_old_fval, gfkp1 = \
                line_search_wolfe12(f, fprime, xk, dk, fpk, old_fval, old_old_fval)
        except _LineSearchError:
            warnflag = 2
            break
        xk1 = xk + step * dk
        sk = xk1 - xk
        xk = xk1
        fpk1 = fprime(xk1)
        yk = fpk1 - fpk
        fpk = fpk1
     # update Hk, DFP formula
        try:
            rho1 = 1.0 / np.dot(sk, yk)
        except ZeroDivisionError:
            rho1 = 1000.0
            print "Divided by Zero,DFP "
        try:
            rho2 = 1.0 / np.dot(yk, np.dot(hk,yk))
        except ZeroDivisionError:
            rho2 = 1000.0
            print "Divided by Zero,DFP "
        if np.isinf(rho1):
            rho1 = 1000.0
            print "nan in DFP"
        if np.isinf(rho2):
            rho2 = 1000.0
            print "nan in DFP"
        hk += (np.outer(sk,sk) * rho1 -
               np.dot(np.outer(np.dot(hk,yk),yk),hk) * rho2)
        iter += 1
    if warnflag == 2:
        print "Line search error in DFP at iteration: ",
        print iter
    else:
        print "DFP Finished"
        print "iterations:"
        print iter
        print 'gk = '
        print fpk
        print 'xk = '
        print xk
#print DFP(watson, np.zeros(4), watson_der, watson_hess)





