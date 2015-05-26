"""
Non-Linear Least Squares Problems
Gauss-Newton Method, LMF, Dogleg
And:
    Large Residual Problems
"""
from scipy.optimize.linesearch import line_search_BFGS, line_search_wolfe1
from scipy.optimize.linesearch import line_search_wolfe2
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

def GN(r, x0, jac, search = True, ave = 1e-8, maxiter = 1000):
    """
    r in vector form
    """
    def f(x):
        return .5 * np.dot(r(x),r(x))
    def fprime(x):
        return np.dot(jac(x),r(x))
    totalfc = 0
    totalgc = 0
    count = 0
    old_fval = f(x0)
    old_old_fval = None
    warnflag = 0
    xk = x0
    jack = jac(xk)
    jack_T = np.transpose(jack)
    dk = LA.solve(np.dot(jack_T,jack), -np.dot(jack_T, jack))
    gk = fprime(xk)
    totalgc += 1
    while (LA.norm(jac(xk)) > ave and count < maxiter):
        if search == True:
            try:
                step, fc, gc, old_fval, old_old_fval, gfkp1 = \
                    line_search_wolfe12(f, fprime, xk, dk, gk, old_fval, old_old_fval)
            except _LineSearchError:
                warnflag = 2
                break
            totalfc += fc
            totalgc += gc
        else:
            step = 1
        count += 1
        xk = xk + step * dk
        jack = jac(xk)
        jack_T = np.transpose(jack)
        dk = LA.solve(np.dot(jack_T,jack), -np.dot(jack_T, jack))
        gk = fprime(xk)
        totalgc += 1
    if count == maxiter:
        warnflag = 3
    if warnflag == 2:
        print "Line search error",
        print " at iteration: ",
        print count
        print_res("Line search error in GN method at iteration:", count, totalfc,
                  totalgc, gk, xk, f(xk))
    elif warnflag == 3:
        print_res("Max number of iterations", count, totalfc,
                  totalgc, gk, xk, f(xk))
    else:
        print_res("Finished", count, totalfc,
                  totalgc, gk, xk, f(xk))













