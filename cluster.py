# Cluster, x = (x1,y1, x2, y2, ..., xn, yn)
# debug exp: sometimes check the math property, not only
# your program
# take care of consistent subscripting!!!
import numpy as np
def u(r):
    return r**(-6) - 2*r**(-3)
def r(x,i,j):
    return (x[2*j] - x[2*i])**2 + (x[2*j+1] - x[2*i+1])**2
def u_der(x,i,j):
    return (-6)*r(x,i,j)**(-7) + 6*r(x,i,j)**(-4)
def cluster(x):
    x = np.asarray(x)
    sum = 0;
    n = len(x) / 2
    for j in range(1, n):
        for i in range(0, j):
            sum +=u(r(x,i,j))
    return sum
def cluster_der(x):
    x = np.asarray(x)
    n = len(x) / 2
    der = np.zeros_like(x)
    #xk: der(2k)
    for k in range(0, n):
        sum = 0
        for j in range(0, k):
            #print k,
            #print j
            sum += (u_der(x,k,j) * 2*(x[2*k] - x[2*j]))
        for j in range(k+1, n):
            sum += (u_der(x,k,j) * 2*(x[2*k] - x[2*j]))
        der[2*k] = sum
        sum = 0
        for j in range(0, k):
            sum += (u_der(x,k,j) * 2*(x[2*k+1] - x[2*j+1]))
        for j in range(k+1, n):
            sum += (u_der(x,k,j) * 2*(x[2*k+1] - x[2*j+1]))
        der[2*k+1] = sum
    return der
