import numpy as np
from watson import watson, watson_der, watson_hess
from propane import propane, propane_der, propane_hess
from cluster import cluster, cluster_der
from dipole import dipole, dipole_der, dipole_hess
import scipy.optimize
import NewtonMethods
from scipy.optimize import rosen, rosen_der, rosen_hess



# -------- Watson ---------
#for i in [3, 6, 9, 12, 15]:
    #print "i = ",
    #print i
    #precision = 1e-10
    #NewtonMethods.SR1(watson, np.zeros(i), watson_der, watson_hess,ave = precision)
    #NewtonMethods.BFGS(watson, np.zeros(i), watson_der, watson_hess, ave = precision)
    #NewtonMethods.DFP(watson, np.ones(i), watson_der, watson_hess)
    #NewtonMethods.stepNewton(watson, np.zeros(i), watson_der, watson_hess,precision)
    #NewtonMethods.adjustedNewton(watson, np.zeros(i), watson_der, watson_hess, precision)

 ##-------- Dipole ---------
#x1 = np.array([.299, .186, -.0273, .0254, -.474, .474, -.892, 0.892])
#x2=np.array([-.3, -.39, .3, -3.44, -1.2, 2.69, 1.59, -1.5])
#x3=np.array([-.041,-.775,.03, -.047, -2.565, 2.565, -.754, .754])
#x1sol = np.array([-6.32e-3, 4.91e-1, -1.99e-3, 9.81e-5, 1.22e-1, -1.00e-1, -4.02, -2.07e-2])
#x2sol = np.array([-3.11e-1, -3.78e-1, 3.28e-1, -3.72e-1, -1.28, 2.49, 1.55, -1.38])
#res = scipy.optimize.minimize(dipole, x3, method='BFGS',jac = dipole_der,
               #options={'disp':True})
#print res.x
#x = x3
#NewtonMethods.stepNewton(dipole, x, dipole_der, dipole_hess)
#NewtonMethods.adjustedNewton(dipole, x, dipole_der, dipole_hess)
#NewtonMethods.SR1(dipole, x, dipole_der, dipole_hess, maxiter = 2000)
#NewtonMethods.BFGS(dipole, x, dipole_der, None, maxiter = 2000)
#NewtonMethods.DFP(dipole, x, dipole_der, dipole_hess)



# --------- Test---------
#x0 = [1.3, 0.7, 0.8, 1.9, 1.2]
#x1 = [1.1, 0.9, 1.1, 1.2, 1.1]
#res = scipy.optimize.minimize(rosen, x0, method='BFGS',jac = rosen_der,
               #options={'disp':True})
#print res.x
#print x0
#NewtonMethods.BFGS(rosen, x0, rosen_der, rosen_hess)
#NewtonMethods.stepNewton(rosen, x0, rosen_der, rosen_hess)
#NewtonMethods.SR1(rosen, x1, rosen_der, rosen_hess)
#NewtonMethods.DFP(rosen, x0, rosen_der, rosen_hess, maxiter = 1500)



# -------- Cluster ---------
def clusterPoints(n):
    m =int(np.sqrt(n / 4))
    pts = np.zeros(2 * n)
    h = 1
    for i in range(0, m * m):
        pts[2*i] = (i % m)*h + h
        pts[2*i+1] = int(i / m)*h + h
    for i in range(m * m, 2 * m * m):
        j = i % m ** 2
        pts[2*i] = -1 * ((j % m)*h + h)
        pts[2*i+1] = int(j / m)*h + h
    for i in range(2 * m**2, 3 * m**2):
        j = i % m ** 2
        pts[2*i] = -1 * ((j % m)*h + h)
        pts[2*i+1] = -1 * (int(j / m)*h + h)
    for i in range(3 * m**2, 4 * m**2):
        j = i % m ** 2
        pts[2*i] = (j % m)*h + h
        pts[2*i+1] = -1 * (int(j / m)*h + h)
    return pts
xclu = clusterPoints(196)
xclu = np.append(xclu,np.array([.5,.5, -.5, .5, -.5,-.5, .5, -.5]))
res = scipy.optimize.minimize(cluster,xclu,method='BFGS', jac = cluster_der,
                              options={'disp' : True})
print res.x
#print xclu
#NewtonMethods.SR1(cluster, xclu, cluster_der, None)
#NewtonMethods.BFGS(cluster, xclu, cluster_der, None)
#NewtonMethods.DFP(cluster, xclu, cluster_der, None)
#NewtonMethods.DFP(cluster, xclu, cluster_der, None)



# ---------Propane---------
#x0 = np.ones(5);
#xsol = np.array([0.31e-2 , 0.345e2, 0.65e-1, .859, .369e-1])
#xs = 6.1 * xsol
#propane(xs)
#propane_der(xs)
#propane_hess(xs)
#for i in [5, 10, 20]:
    #xs = i * xsol
    #NewtonMethods.BFGS_Modified(propane, xs, propane_der,propane_hess)
    #NewtonMethods.adjustedNewton(propane, xs, propane_der, propane_hess)
    #NewtonMethods.stepNewton(propane,xs,propane_der, propane_hess)
    #NewtonMethods.DFP(propane, xs, propane_der, propane_hess, maxiter = 1000)
    #NewtonMethods.BFGS(propane,xs,propane_der, propane_hess)
    #NewtonMethods.SR1(propane, xs, propane_der, propane_hess)
##print SR1(propane, x0, propane_der, propane_hess, 1e-14)
#res = scipy.optimize.minimize(propane, xs, method='BFGS',jac = propane_der,
                              #options={'disp':True})
#print (res.x)

