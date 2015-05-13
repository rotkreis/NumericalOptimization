from scipy.optimize.linesearch import line_search_BFGS, line_search_wolfe1
from scipy.optimize.linesearch import line_search_wolfe2
from scipy.optimize import rosen, rosen_der, rosen_hess
import numpy as np
from scipy import linalg as LA
import scipy.optimize
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
def cgsolve(f, x0, fprime, method, ave = 1e-8, maxiter = 2000):
    gk = fprime(x0)
    dk = -gk
    xk = x0
    totalfc = 0
    totalgc = 0
    count = 0
    old_fval = f(x0)
    old_old_fval = None
    warnflag = 0
    while (LA.norm(gk) > ave and count < maxiter):
        try:
            step, fc, gc, old_fval, old_old_fval, gfkp1 = \
                line_search_wolfe12(f, fprime, xk, dk, gk, old_fval, old_old_fval)
        except _LineSearchError:
            warnflag = 2
            break
        totalfc += fc
        totalgc += gc
        xk1 = xk + step * dk
        gk1 = fprime(xk1)
        gk = fprime(xk)
        beta = 0
        if method =="FR":
            beta = np.dot(gk1,gk1) / np.dot(gk,gk)
        elif method =="PRP":
            beta = np.dot(gk1, gk1-gk) / np.dot(gk,gk)
        elif method =="PRP+":
            beta = np.abs(np.dot(gk1, gk1-gk) / np.dot(gk,gk))
        elif method == "CD":
            beta = -np.dot(gk1,gk1) / np.dot(dk,gk)
        else:
            # "DY"
            beta = np.dot(gk1,gk1) / np.dot(dk, gk1-gk)
        dk = - gk1 + beta * dk
        xk = xk1
        gk = gk1
        count += 1
    if count == maxiter:
        warnflag = 3
    if warnflag == 2:
        print "Line search error in ",
        print method,
        print " at iteration: ",
        print count
        print_res(method, count, totalfc,totalgc, gk, xk, f(xk))
    elif warnflag == 3:
        print_res(method+" reached max number of iterations",count,totalfc,totalgc,gk,xk,f(xk))
    else:
        print_res(method+"Finished",count,totalfc,totalgc,gk,xk,f(xk))

