from scipy.optimize.linesearch import line_search_wolfe2 as line_search
import numpy as np
from scipy import linalg as LA
import scipy.optimize
from watson import watson, watson_der, watson_hess

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
        xk += dk * line_search(f, fprime, xk, dk, fpk)[0]
        iter += 1
        fpk = fprime(xk)
    print "iterations:"
    print iter
    print LA.norm(fpk)
    return f(xk)

#print stepNewton(watson, np.zeros(6), watson_der, watson_hess)

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
        xk += dk * line_search(f, fprime, xk, dk, fpk)[0]
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
        step = line_search(f, fprime, xk, dk, fpk)[0]
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
    hk = LA.inv(fhess(x0))
    while (LA.norm(fpk) > ave and iter < maxiter):
        dk = - np.dot(hk,fpk)
        step = line_search(f, fprime, xk, dk, fpk)[0]
        if step == None:
            print "line search error"
        sk = dk * step
        xk += sk
        yk = fprime(xk) - fpk
        fpk += yk
        # update Hk, BFGS formula
        hk += ((1 + np.dot(yk, np.dot(hk,yk))/np.dot(yk,sk)) * np.outer(sk,sk) / np.dot(yk,sk)
               -(np.dot(np.outer(sk,yk),hk) + np.dot(hk,np.outer(yk,sk))) / np.dot(yk,sk))
        iter += 1
    print "iterations:"
    print iter
    print f(xk)
    return f(xk)
print BFGS(watson, np.zeros(9), watson_der, watson_hess)


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
        step = line_search(f, fprime, xk, dk, fpk)[0]
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
    print f(xk)
    return f(xk)
#print DFP(watson, np.zeros(9), watson_der, watson_hess, 1e-5)














