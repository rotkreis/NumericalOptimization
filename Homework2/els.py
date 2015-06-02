import numpy as np
golden = (-1 + np.sqrt(5))/2
def interval(f,d,x0, a0 = 1, y0 = 1, t = golden+1):
    count = 0
    ai = a0
    yi= y0
    while True:
        a1 = ai+ yi
        if a1 <= 0 or f(x0 + a1*d) >= f(x0 + ai*d):
            if a1 <= 0:
                a1 = 0
            if count == 0:
                yi= -yi
                a = a1
            else:
                alpha = np.min([a,a1])
                beta = np.max([a,a1])
                return alpha, beta
        yi = t * yi
        a = ai
        ai = a1
        count += 1

def f(x):
    return np.cos(x)
    #return (x-1)**2
    #if x <= 2:
        #return (x-2)**2
    #elif x >= 2:
        #return (x-4)**2

[a,b] = interval(f, 1, 0)
print [a,b]

def goldensection(f, d, x0, t = golden, ave = 1e-6, maxiter = 200):
    [ai,bi] = interval(f,d,x0)
    count = 0
    while bi - ai > ave and count <= maxiter:
        al = ai + (1-t)*(bi-ai)
        ar = ai + t*(bi-ai)
        if f(x0+al*d) < f(x0+ar*d):
            bi = ar
        else:
            ai = al
        count += 1
    return (ai+bi)/2

print goldensection(f, 1, 0)

