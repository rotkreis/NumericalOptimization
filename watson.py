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

def der_r(i,x,k):
    if k == 1:
        return -2 * np.dot(x, t(i) ** (np.arange(1, len(x) + 1) - 1))
    elif k == 2:
        return 1 - 2 * np.dot(x, t(i) ** (np.arange(1, len(x) + 1) - 1)) * t(i)
    else:
        return ((k - 1) * t(i) ** (k - 2) -
               2 * np.dot(x, t(i) ** (np.arange(1, len(x) + 1) - 1)) * t(i) ** (k - 1))

def der_watson(x):
    der = np.zeros_like(x)
    # der[0]
    sum = 0
    for i in range(1, 29 + 1):
        sum += 2 * r(i,x) * der_r(i, x, 1)
    sum += 2 * x[0] + 2 * (x[1] - x[0] ** 2 - 1) * (-2 * x[0])
    der[0] = sum
    # der[1]
    sum = 0
    for i in range(1, 29 + 1):
        sum += 2 * r(i,x) * der_r(i, x, 2)
    sum += 2 * (x[1] - x[0] ** 2 - 1)
    der[1] = sum
    # der[2] ~ der[len(der) - 1]
    for k in range(2, len(x)):
        sum = 0
        for i in range(1, 29 + 1):
            sum += 2 * r(i,x) * der_r(i, x, k + 1)
        der[k] = sum
    return der
test_array = np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7],float)
print watson(np.arange(1,7))
print scipy.optimize.check_grad(watson,der_watson, test_array)
print der_watson(test_array)
print scipy.optimize.approx_fprime(test_array, watson,0.000000001)
