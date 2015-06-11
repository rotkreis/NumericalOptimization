from scipy.optimize.linesearch import line_search_BFGS,line_search_armijo, line_search_wolfe1
from scipy.optimize.linesearch import line_search_wolfe2
from scipy.optimize import rosen, rosen_der, rosen_hess
import numpy as np
from scipy import linalg as LA
import scipy.optimize
#from watson import watson, watson_der, watson_hess
#from propane import propane, propane_der, propane_hess
#from cluster import cluster, cluster_der
#hah
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
    totalfc = 0
    totalgc = 0
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
        totalfc += fc
        totalgc += gc
        xk = xk + dk * step
        fpk = fprime(xk)
        iter += 1
    if iter == maxiter:
        warnflag = 3
    if warnflag == 2:
        print "Line search error in stepNewton at iteration: ",
        print iter
        print_res("",iter,totalfc,totalgc, fpk, xk, f(xk))
    elif warnflag == 3:
        print_res("Reached max number of iterations", iter,totalfc, totalgc,  fpk, xk, f(xk))
    else:
        print_res("StepNewton Finished", iter,totalfc, totalgc,  fpk, xk, f(xk))

#stepNewton(watson, np.zeros(9), watson_der, watson_hess, 1e-8)

def adjustedNewton(f, x0, fprime, fhess, u = 1e-4, ave = 1e-6, maxiter = 2000):
    """
    my simple minimization algorithm, testing
    using LM - newton Method
    """
    totalfc = 0
    totalgc = 0
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
                hess = hess + u * np.eye(len(x0),len(x0))
                u *= 2
        dk = LA.solve(hess, fpk * -1)
        try:
            step, fc, gc, old_fval, old_old_fval, gfkp1 = \
                line_search_wolfe12(f, fprime, xk, dk, fpk, old_fval, old_old_fval)
        except _LineSearchError:
            warnflag = 2
            break
        totalfc += fc
        totalgc += gc
        xk = xk + dk * step
        iter += 1
        fpk = fprime(xk)
    if iter == maxiter:
        warnflag = 3
    if warnflag == 2:
        print "Line search error in adjustedNewton at iteration: ",
        print iter
        print_res("",iter,totalfc,totalgc, fpk, xk, f(xk))
    elif warnflag == 3:
        print_res("Reached max number of iterations", iter,totalfc, totalgc,  fpk, xk, f(xk))
    else:
        print_res("Adjusted Newton Finished", iter,totalfc, totalgc,  fpk, xk, f(xk))
#adjustedNewton(watson, np.zeros(9), watson_der, watson_hess, 1e-8)

def SR1(f, x0, fprime, fhess, ave = 1e-6, maxiter = 1000):
    """
    Quasi-Newton: SR1
    fhess: inital hessian, converted to H0 for further compuatitions
    possible error: the update of sk, yk, hk
    """
    iter = 0
    totalfc = 0
    totalgc = 0
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
        totalfc += fc
        totalgc += gc
        sk = dk * step
        xk = xk + sk
        fpk1 = fprime(xk)
        yk = fpk1 - fpk
        fpk = fpk1
        temp = sk - np.dot(hk,yk)
        try:
            rhok = 1.0 / np.dot(temp, yk)
        except ZeroDivisionError:
            rhok = 1000.0
            if warnflag == 0:
                warnflag = 1
                print "Divided by Zero, SR1"
        if np.isinf(rhok):
            rhok = 1000.0
            if warnflag == 0:
                warnflag = 1
                print "Inf in SR1"
        hk = hk + (np.outer(temp, temp) * rhok)
        iter += 1

    if iter == maxiter:
        warnflag = 3
    if warnflag == 2:
        print "Line search error in SR1 at iteration: ",
        print iter
        print_res("",iter,totalfc,totalgc, fpk, xk, f(xk))
    elif warnflag == 1:
        print_res("SR1 questionable", iter,totalfc, totalgc,  fpk, xk, f(xk))
    elif warnflag == 3:
        print_res("Reached max number of iterations", iter,totalfc, totalgc,  fpk, xk, f(xk))
    else:
        print_res("SR1 Finished",iter,totalfc, totalgc, fpk,xk, f(xk))

#SR1(watson, np.ones(9), watson_der, watson_hess)
def BFGS(f, x0, fprime, fhess = None, ave = 1e-6, maxiter = 1000,
         diag = False, disp = False, wolfe2 = .25):
    """
    Quasi-Newton: BFGS
    fhess: inital hessian, converted to H0 for further compuatitions
    possible error: the update of sk, yk, hk
    """
    iter = 0
    totalfc = 0
    totalgc = 0
    fpk = fprime(x0)
    xk = x0
    old_fval = f(x0)
    old_old_fval = None
    I = np.eye(len(x0))
    # Holy Shit! Initially inv(hessian!) now ok ! hahahhahaha
    hk = I
    warnflag = 0
    while (LA.norm(fpk) > ave and iter < maxiter):
        dk = - np.dot(hk,fpk)
        try:
            step, fc, gc, old_fval, old_old_fval, gfkp1 = \
                line_search_wolfe12(f, fprime, xk, dk, fpk, old_fval, old_old_fval, c2 = wolfe2)
        except _LineSearchError:
            warnflag = 2
            break
        totalfc += fc
        totalgc += gc
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

    msg = " "
    if disp == True:
            if iter == maxiter:
                warnflag = 3
            if warnflag == 2:
                msg = "Line search error in BFGS "
                print "Line search error in BFGS at iteration: ",
                print iter
                print_res("",iter,totalfc,totalgc, fpk, xk, f(xk))
            elif warnflag == 3:
                msg = "BFGS max number of iterations"
                print_res(msg, iter,totalfc, totalgc,  fpk, xk, f(xk))
            else:
                print_res("BFGS Finished",iter,totalfc, totalgc, fpk,xk, f(xk))

    return iter,totalfc, totalgc, fpk,xk, f(xk), warnflag, msg
    #return xk,warnflag, msg, iter, totalfc, totalgc, f(xk), fprime(xk)
#BFGS(watson, np.zeros(9), watson_der, watson_hess)

def DFP(f, x0, fprime, fhess, ave = 1e-6, maxiter = 1000):
    """
    Quasi-Newton: DFP
    fhess: inital hessian, converted to H0 for further compuatitions
    possible error: the update of sk, yk, hk
    """
    iter = 0
    totalfc = 0
    totalgc = 0
    fpk = fprime(x0)
    xk = x0
    old_fval = f(x0)
    old_old_fval = None
    I = np.eye(len(x0), dtype=float)
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
        totalfc += fc
        totalgc += gc
        xk1 = xk + step * dk
        sk = xk1 - xk
        xk = xk1
        fpk1 = fprime(xk1)
        yk = fpk1 - fpk
        fpk = fpk1
        try:
            rho1 = 1.0 / np.dot(sk, yk)
        except ZeroDivisionError:
            rho1 = 1000.0
            if warnflag == 0:
                print "Divided by Zero,DFP at iteration ",
                print iter
                warnflag = 1
        try:
            rho2 = 1.0 / np.dot(yk, np.dot(hk,yk))
        except ZeroDivisionError:
            rho2 = 1000.0
            if warnflag == 0:
                print "Divided by Zero,DFP at iteration ",
                print iter
                warnflag = 1
        if np.isinf(rho1):
            rho1 = 1000.0
            if warnflag == 0:
                print "Inf in DFP at iteration ",
                print iter
                warnflag = 1
        if np.isinf(rho2):
            rho2 = 1000.0
            if warnflag == 0:
                print "Inf in DFP at iteration ",
                print iter
                warnflag = 1
        hk = hk + (np.outer(sk,sk) * rho1 -
                   np.dot(np.outer(np.dot(hk,yk),yk),hk) * rho2)
        iter += 1

    if iter == maxiter:
        warnflag = 3
    if warnflag == 2:
        print "Line search error in DFP at iteration: ",
        print iter
        print_res("",iter,totalfc,totalgc, fpk, xk, f(xk))
    elif warnflag == 1:
        print_res("DFP questionable",iter,totalfc,totalgc, fpk, xk, f(xk))
    elif warnflag == 3:
        print_res("Reached max number of iterations", iter,totalfc, totalgc,  fpk, xk, f(xk))
    else:
        print_res("DFP Finished",iter,totalfc, totalgc, fpk,xk, f(xk))

#DFP(watson, np.ones(4), watson_der, watson_hess)
x0 = [1.3, 0.7, 0.8, 1.9, 1.2]
#res = scipy.optimize.minimize(rosen, x0, method='BFGS',jac = rosen_der,
               #options={'disp':True})
#print res.x
#DFP(rosen, x0, rosen_der, rosen_hess)

