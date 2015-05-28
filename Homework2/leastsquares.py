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
def Dogleg(r, x0, jac, delta0 = 1, ave = 1e-8, maxiter = 1000):
    """
    Trust Region Method
    """
    def f(x):
        return .5 * np.dot(r(x),r(x))
    def fprime(x):
        return np.dot(np.transpose(jac(x)),r(x))
    totalfc = 0
    totalgc = 0
    count = 0
    warnflag = 0
    xk = x0
    deltak = delta0
    dim = x0.size
    while (LA.norm(jac(xk)) > ave and count < maxiter):
        jack = jac(xk)
        jack_T = np.transpose(jack)
        rk = r(xk)
        dGNk = LA.solve(np.dot(jack_T,jack), -np.dot(jack_T, rk))
        dSDk = - np.dot(np.transpose(jack), rk)
        alphak = LA.norm(dSDk)**2 / LA.norm(np.dot(jack,dSDk))**2
        if LA.norm(dGNk) <= deltak:
            dk = dGNk
        elif alphak * LA.norm(dSDk) >= deltak:
            dk = deltak * dSDk / LA.norm(dSDk)
        else:
            c1 = LA.norm(dGNk - alphak * dSDk)**2
            c2 = 2 * alphak * np.dot(dSDk, dGNk - alphak * dSDk)
            c3 = alphak ** 2 * LA.norm(dSDk)**2 - deltak**2
            [b1,b2] = np.roots([c1,c2,c3])
            if b1 > 1e-10:
                beta = b1
            dk = (1-beta) * alphak * dSDk + beta * dGNk
        xk1 = xk + dk
        dfk = f(xk) - f(xk1)
        def qk(d):
            return .5 * LA.norm(np.dot(jack, dk) + rk)**2
        dqk = qk(np.zeros(dim)) - qk(dk)
        gammak = dfk / dqk
        if gammak > 0.75 and LA.norm(dk) == deltak:
            deltak = 2 * deltak
        elif gammak < 0.25:
            deltak = 0.25 * deltak
        if gammak > 0:
            xk = xk1
            count += 1
    print_res("Dogleg Finished", count, totalfc, totalgc, fprime(xk), xk, f(xk))

def LMF(r, x0, jac, v0 = 1e-2, ave = 1e-8, maxiter = 1000):
    def f(x):
        return .5 * np.dot(r(x),r(x))
    def fprime(x):
        return np.dot(np.transpose(jac(x)),r(x))
    totalfc = 0
    totalgc = 0
    count = 0
    warnflag = 0
    xk = x0
    vk = v0
    dim = x0.size
    while (LA.norm(jac(xk)) > ave and count < maxiter):
        jack = jac(xk)
        jack_T = np.transpose(jack)
        dk = LA.solve(np.dot(jack_T,jack) + vk*np.eye(dim), -np.dot(jack_T, r(xk)))
        gk = fprime(xk)
        totalgc += 1
        xk1 = xk + dk
        dfk = f(xk) - f(xk1)
        dqk = .5 * np.dot(dk, vk*dk - gk)
        gammak = dfk / dqk
        if gammak < 0.25:
            vk = 4*vk
        elif gammak > 0.75:
            vk = vk / 2.0
        if gammak > 0:
            xk = xk + dk
            count += 1
    print_res("LMF Finished", count, totalfc, totalgc, gk, xk, f(xk))


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
    dk = LA.solve(np.dot(jack_T,jack), -np.dot(jack_T, r(xk)))
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
        dk = LA.solve(np.dot(jack_T,jack), -np.dot(jack_T, r(xk)))
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













