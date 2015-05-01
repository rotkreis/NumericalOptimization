#from scipy.optimize.linesearch import line_search_wolfe2 as line_search
from cluster import cluster, cluster_der
import numpy as np
import scipy.optimize
from dipole import dipole, dipole_der, dipole_hess
#from scipy import linalg as LA
#import scipy.optimize

#def propane(x1,x2,x3,x4,x5,p = 40, R = 10, r5,r6,r7,r8,r9,r10):
    ##define constant parameters k_n
    #f = np.zeros(5)
    #f[0] = x1 * x2 + x1 - 3 * x5
    #f[1] = (2*x1*x2 + x1+ 2*r10*x2**2 + x2*x3**2 + r7*x2*x3+
            #r9*x2*x4 + r8*x2 - R*x5)
    #f[2] = 2*x2*x3**2 + r7*x2*x3 + 2*r5*x3**2 + r6*x3 - 8*x5
    #f[3] = r9*x2*x4 + 2*x4**2 - 4*R*x5
    #f[4] = (x1*x2 + x1 + r10*x2**2 + x2*x3**2 + r7*x2*x3 + r9*x2*x4 + r8*x2 +
            #r5*x3**2 + r6*x3 + x4**2 - 1)
    #return LA.norm(f,2)
#def propane_der(x1,x2,x3,x4,x5,R,r5,r6,r7,r8,r9,r10):
#print scipy.optimize.check_grad(cluster, cluster_der,x0)
#print scipy.optimize.approx_fprime(x0,cluster,1e-8)

x1 = np.array([.299, .186, -.0273, .0254, -.474, .474, -.892, 0.892])
x1sol = np.array([-6.32e-3, 4.91e-1, -1.99e-3, 9.81e-5, 1.22e-1, -1.00e-1, -4.02, -2.07e-2])
x2sol = np.array([-3.11e-1, -3.78e-1, 3.28e-1, -3.72e-1, -1.28, 2.49, 1.55, -1.38])
print scipy.optimize.check_grad(dipole, dipole_der, x1sol)
print dipole_der(x1sol)
res = scipy.optimize.minimize(dipole, x1, method='BFGS',jac = dipole_der,
               options={'disp':True})
print res.x

