import numpy as np
import scipy.optimize
def t(i):
    return i / 29.0

def r(i,x):
    sum = -1
    sum += np.dot(x * (np.arange(1, len(x) + 1) - 1),t(i) ** (np.arange(1,len(x) + 1) - 2))
    sum -= np.dot(x, t(i) ** (np.arange(1,len(x)+1) - 1)) ** 2
    return sum

def watson(x):
    sum = 0
    for i in range(1,29 + 1):
        sum += r(i,x) ** 2
    sum += x[0] ** 2
    sum += (x[1] - x[0] ** 2 - 1) ** 2
    return sum

def r_der(i,x,k):
    if k == 1:
        return -2 * np.dot(x, t(i) ** (np.arange(1, len(x) + 1) - 1))
    elif k == 2:
        return 1 - 2 * np.dot(x, t(i) ** (np.arange(1, len(x) + 1) - 1)) * t(i)
    else:
        return ((k - 1) * t(i) ** (k - 2) -
               2 * np.dot(x, t(i) ** (np.arange(1, len(x) + 1) - 1)) * t(i) ** (k - 1))

def watson_der(x):
    der = np.zeros_like(x)
    # der[0]
    sum = 0
    for i in range(1, 29 + 1):
        sum += 2 * r(i,x) * r_der(i, x, 1)
    sum += 2 * x[0] + 2 * (x[1] - x[0] ** 2 - 1) * (-2 * x[0])
    der[0] = sum
    # der[1]
    sum = 0
    for i in range(1, 29 + 1):
        sum += 2 * r(i,x) * r_der(i, x, 2)
    sum += 2 * (x[1] - x[0] ** 2 - 1)
    der[1] = sum
    # der[2] ~ der[len(der) - 1]
    for k in range(2, len(x)):
        sum = 0
        for i in range(1, 29 + 1):
            sum += 2 * r(i,x) * r_der(i, x, k + 1)
        der[k] = sum
    return der
test_array = np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7],float)
#print watson(np.arange(1,7))
#print scipy.optimize.check_grad(watson,watson_der, test_array)
#print watson_der(test_array)
#print scipy.optimize.approx_fprime(test_array, watson,0.000000001)

def r_der(i,x,k,l):
    return - 2 * t(i) ** ( k + l - 2)

def watson_hess(x):
    hess = np.zeros((len(x),len(x)))
#    # hess[0,0]
#    sum = 0
#    for i in range(1, 29 + 1):
#        sum += 2 * r_der(i,x,1) ** 2 + 2 * r(i) * r_der(i,x,1,1)
#    sum += 2 + 12 * x[0] ** 2 - 4 * x[1] + 4
#    hess[0,0] = sum
#    # hess[0,1], hess[1,0]
#    sum = 0
#    for i in range(1, 29 + 1):
#        sum += 2 * r_der(i,x,2) * r_der(i,x,1) + 2 * r(i) * r_der(i,x,1,2)
#    sum += -4 * x[0]
#    hess[0,1] = hess[1,0] = sum
#    # hess[1,1]
#    sum = 0
#    for i in range(1, 29 + 1):
#        sum += 2 * r_der(i,x,2) ** 2 + 2 * r(i) * r_der(i,x,2,2)
#    sum += 2
#    hess[1,1] = sum
    for k in range(0,len(x)):
        for l in range(0, len(x)):
            sum = 0
            for i in range(1,29 + 1):
                sum += (2 * r_der(i,x,k) * r_der(i, x, l) +
                        2 * r(i) * r_der(i,x,k,l))
    hess[0,0] += 2 + 12 * x[0] ** 2 - 4 * x[1] + 4
    hess[0,1] += -4 * x[0]
    hess[1,0] = hess[0,1]
    hess[1,1] += 2
    return hess






















