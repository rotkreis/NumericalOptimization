from scipy.optimize.linesearch import line_search_BFGS,line_search_armijo, line_search_wolfe1
from scipy.optimize.linesearch import line_search_wolfe2
from scipy.optimize import rosen, rosen_der, rosen_hess
import numpy as np
from scipy import linalg as LA
import scipy.optimize
from watson import watson, watson_der, watson_hess
from propane import propane, propane_der, propane_hess
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
        print "Wolfe12 Error"
    return ret

def stepNewton(f, x0, fprime, fhess, ave = 1e-5, maxiter = 2000):
    """
    my simple minimization algorithm, testing
    using newton Method
    """
    iter = 0
    fpk = fprime(x0)
    xk = x0
    print ave
    while (LA.norm(fpk) > ave and iter < maxiter):
        hess = fhess(xk)
        dk = LA.solve(hess, fpk * -1)
        xk += dk * line_search_wolfe12(f, fprime, xk, dk, fpk)[0]
        iter += 1
        fpk = fprime(xk)
    print "iterations:"
    print iter
    print LA.norm(fpk)
    return f(xk)

xs = np.array([0.31e-2 , 0.345e2, 0.65e-1, .859, .369e-1])
xs = 10.01 * xs
#print propane(xs)
#print propane_der(xs)
#print propane_hess(xs)
#print stepNewton(watson, np.zeros(6), watson_der, watson_hess)
#print stepNewton(propane,xs,propane_der, propane_hess)

def adjustedNewton(f, x0, fprime, fhess, u = 1e-4, ave = 1e-6, maxiter = 2000):
    """
    my simple minimization algorithm, testing
    using LM - newton Method
    """
    iter = 0
    fpk = fprime(x0)
    xk = x0
    while (LA.norm(fpk) > ave and iter < maxiter):
        hess = fhess(xk)
        l = LA.cholesky(hess)
        while any(np.diag(l) < 0 ):
            hess += u * np.eye(len(x0),len(x0))
            u *= 2
        dk = LA.solve(hess, fpk * -1)
        xk += dk * line_search_wolfe12(f, fprime, xk, dk, fpk)[0]
        iter += 1
        fpk = fprime(xk)
    print "iterations:"
    print iter
    return f(xk)
#print adjustedNewton(watson, np.zeros(9), watson_der, watson_hess, 1e-8)

def SR1(f, x0, fprime, fhess, ave = 1e-6, maxiter = 1000):
    """
    Quasi-Newton: SR1
    fhess: inital hessian, converted to H0 for further compuatitions
    possible error: the update of sk, yk, hk
    """
    iter = 0
    fpk = fprime(x0)
    xk = x0
    hk = LA.inv(fhess(x0))
    while (LA.norm(fpk) > ave and iter < maxiter):
        dk = - np.dot(hk,fpk)
        step = line_search_wolfe12(f, fprime, xk, dk, fpk)[0]
        if step == None:
            print "line search error"
        sk = dk * step
        xk += sk
        yk = fprime(xk) - fpk
        fpk += yk
        temp = sk - np.dot(hk,yk)
        hk += (np.outer(temp, temp) /
               np.dot(temp, yk))
        iter += 1
    print "iterations:"
    print iter
    print f(xk)
    return f(xk)

#print SR1(watson, np.zeros(9), watson_der, watson_hess)

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
    old_old_fval = None
    I = np.eye(len(x0), dtype=int)
    # Holy Shit! Initially inv(hessian!) now ok ! hahahhahaha
    hk = I
    while (LA.norm(fpk) > ave and iter < maxiter):
        dk = - np.dot(hk,fpk)
        step = line_search_wolfe12(f, fprime, xk, dk, fpk, old_fval, old_old_fval)[0]
        #step = line_search_BFGS(f, xk, dk, fpk, k)[0]
        if step == None:
            print "line search error"
        #sk = dk * step
        #xk += sk
        #yk = fprime(xk) - fpk
        #fpk += yk
        xk1 = xk + step * dk
        sk = xk1 - xk
        xk = xk1
        fpk1 = fprime(xk1)
        yk = fpk1 - fpk
        fpk = fpk1
        # update Hk, BFGS formula
        #hk += ((1 + np.dot(yk, np.dot(hk,yk))/np.dot(yk,sk)) * np.outer(sk,sk) / np.dot(yk,sk)
               #-(np.dot(np.outer(sk,yk),hk) + np.dot(hk,np.outer(yk,sk))) / np.dot(yk,sk))
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
    print "iterations:"
    print iter
    print xk
    return f(xk)
#print BFGS(watson, np.zeros(9), watson_der, watson_hess)
#xs = [0.01,100,0.1,0.1,0.1]

def DFP(f, x0, fprime, fhess, ave = 1e-6, maxiter = 1000):
    """
    Quasi-Newton: DFP
    fhess: inital hessian, converted to H0 for further compuatitions
    possible error: the update of sk, yk, hk
    """
    iter = 0
    fpk = fprime(x0)
    xk = x0
    hk = LA.inv(fhess(x0))
    while (LA.norm(fpk) > ave and iter < maxiter):
        dk = - np.dot(hk,fpk)
        step = line_search_wolfe12(f, fprime, xk, dk, fpk)[0]
        if step == None:
            print "line search error"
        sk = dk * step
        xk += sk
        yk = fprime(xk) - fpk
        fpk += yk
        # update Hk, DFP formula
        hk += (np.outer(sk,sk)/np.dot(sk,yk) -
               np.dot(np.outer(np.dot(hk,yk), yk), hk) /
               np.dot(yk, np.dot(hk,yk)))
        iter += 1
    print "iterations:"
    print iter
    return f(xk)
#print DFP(watson, np.zeros(9), watson_der, watson_hess, 1e-5)
#print DFP(propane, xs, propane_der, propane_hess)
#print adjustedNewton(propane,xs,propane_der, propane_hess)
print BFGS(propane,xs,propane_der, propane_hess)
#print SR1(propane, xs, propane_der, propane_hess, 1e-14)
#res = scipy.optimize.minimize(propane, xs, method='BFGS',jac = propane_der,
                              #options={'disp':True})
#print (res.x)
#x0 = [1.3, 0.7, 0.8, 1.9, 1.2]
#res = scipy.optimize.minimize(rosen, x0, method='BFGS',jac = rosen_der,
               #options={'disp':True})
#print res.x
#print x0
#print BFGS(rosen, x0, rosen_der, rosen_hess)












