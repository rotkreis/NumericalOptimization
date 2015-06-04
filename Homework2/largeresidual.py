"""
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

def DGW(r, x0, jac, h = None, search = True, ave = 1e-8, maxiter = 1000,
        diag = False, method = None):
    """
    r in vector form
    """
    def f(x):
        return .5 * np.dot(r(x),r(x))
    def fprime(x):
        return np.dot(np.transpose(jac(x)),r(x))
    totalfc = 0
    totalgc = 0
    count = 0
    old_fval = f(x0)
    old_old_fval = None
    warnflag = 0
    xk = x0
    jk = jac(xk)
    jk_T = np.transpose(jk)
    rk = r(xk)
    dim = len(x0)
    if h == None:
        bk = np.eye(dim) * 1e-2
    else:
        bk = h(x0)
    # question, how to choose the initial value of bk?
    # and whether add it here or not?
    dk = LA.solve(np.dot(jk_T,jk)+bk , -np.dot(jk_T, rk))
    gk = fprime(xk)
    totalgc += 1
    errcount = 0
    while (LA.norm(fprime(xk)) > ave and count < maxiter):
        if search == True:
            if np.dot(dk,gk) > 0:
                dk = -dk
            try:
                step, fc, gc, old_fval, old_old_fval, gfkp1 = \
                    line_search_wolfe12(f, fprime, xk, dk, gk, old_fval, old_old_fval)
            except _LineSearchError:
                print "descend direction:",
                print np.dot(dk,gk)
                warnflag = 2
                break
                #step = 1e-1
                #errcount = errcount + 1
                #if errcount >= 5:
                    #break
            totalfc += fc
            totalgc += gc
        else:
            step = 1
        xk1 = xk + step * dk
        jk1 = jac(xk1)
        jk1_T = np.transpose(jk1)
        rk1 = r(xk1)
        if diag == True:
            print count
            print rk
            print jk
            print xk
            print totalfc
            print f(xk)

        # modify bk
        sk = step * dk
        yk = np.dot(jk1_T, rk1) - np.dot(jk_T, rk)
        ykhat = np.dot(jk1_T,rk1) - np.dot(jk_T, rk1)
        if method == None: # DGW method
            A = np.outer(ykhat - np.dot(bk,sk), yk) + np.outer(yk, ykhat - np.dot(bk,sk))
            B = np.dot(ykhat - np.dot(bk,sk),sk) * np.outer(yk,yk)
            # gamma: adjust convergence
            gamma0 = np.dot(sk,ykhat)/np.dot(sk, np.dot(bk,sk))
            gamma = np.min([1, np.abs(gamma0)])
            bk = gamma * bk
            bk = bk + A / np.dot(yk,sk) - B / np.dot(yk,sk)**2
        elif method == "Biggs":
            gamma = np.dot(rk1,rk1) / np.dot(rk,rk)
            bk = gamma * bk
            A = ykhat - np.dot(bk,sk)
            bk = bk + np.outer(A,A) / np.dot(A,sk)
        elif method == "BFGS":
            bk = bk + np.outer(ykhat,ykhat)/np.dot(ykhat,sk)
            A = np.dot(bk,sk)
            bk = bk - np.outer(A,A) / np.dot(sk,A)

        dk = LA.solve(np.dot(jk1_T, jk1) + bk, -np.dot(jk1_T, rk1))
        totalfc += 1
        gk = fprime(xk1)
        if LA.norm(f(xk1)-f(xk)) < ave:
            break
        xk = xk1
        jk = jk1
        jk_T = jk1_T
        rk = rk1
        totalgc += 1
        count += 1

    if count == maxiter:
        warnflag = 3
    if warnflag == 2:
        print "Line search error",
        print " at iteration: ",
        print count
        print_res("Line search error in DGW method at iteration:", count, totalfc,
                  totalgc, gk, xk, f(xk))
    elif warnflag == 3:
        print_res("DGW Max number of iterations", count, totalfc,
                  totalgc, gk, xk, f(xk))
    else:
        print_res("DGW Finished", count, totalfc,
                  totalgc, gk, xk, f(xk))
    #print LA.norm(fprime(xk))
