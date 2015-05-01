import numpy as np
from watson import watson, watson_der, watson_hess
from propane import propane, propane_der, propane_hess
from cluster import cluster, cluster_der
from dipole import dipole, dipole_der, dipole_hess
import scipy.optimize
import NewtonMethods
from scipy.optimize import rosen, rosen_der, rosen_hess

# -------- Dipole ---------
x1 = np.array([.299, .186, -.0273, .0254, -.474, .474, -.892, 0.892])
x1sol = np.array([-6.32e-3, 4.91e-1, -1.99e-3, 9.81e-5, 1.22e-1, -1.00e-1, -4.02, -2.07e-2])
x2sol = np.array([-3.11e-1, -3.78e-1, 3.28e-1, -3.72e-1, -1.28, 2.49, 1.55, -1.38])
#res = scipy.optimize.minimize(dipole, x1, method='BFGS',jac = dipole_der,
               #options={'disp':True})
#print res.x
#NewtonMethods.adjustedNewton(dipole, x1, dipole_der, dipole_hess)
#NewtonMethods.stepNewton(dipole, x1, dipole_der, dipole_hess)
#NewtonMethods.DFP(dipole, x1, dipole_der, dipole_hess)
#NewtonMethods.SR1(dipole, x1, dipole_der, dipole_hess, maxiter = 2000)
#NewtonMethods.BFGS(dipole, x1, dipole_der, None)



# --------- Test---------
x0 = [1.3, 0.7, 0.8, 1.9, 1.2]
#res = scipy.optimize.minimize(rosen, x0, method='BFGS',jac = rosen_der,
               #options={'disp':True})
#print res.x
#print x0
#NewtonMethods.BFGS(rosen, x0, rosen_der, rosen_hess)
#NewtonMethods.stepNewton(rosen, x0, rosen_der, rosen_hess)
x1 = [1.1, 0.9, 1.1, 1.2, 1.3]
NewtonMethods.SR1(rosen, x1, rosen_der, rosen_hess)
#NewtonMethods.DFP(rosen, x0, rosen_der, rosen_hess, maxiter = 1500)



# -------- Cluster ---------
#xclu=np.array([-5.0,0.0, 0.0,1.0, 1.0,2.0])
#print NewtonMethods.BFGS(cluster, xclu, cluster_der, None)
#print NewtonMethods.DFP(cluster, xclu, cluster_der, None)



# ---------Propane---------
#x0 = np.ones(5);
#xs = np.array([0.31e-2 , 0.345e2, 0.65e-1, .859, .369e-1])
#xs = 6.1 * xs
#print propane(xs)
#print propane_der(xs)
#print propane_hess(xs)
#print stepNewton(watson, np.zeros(6), watson_der, watson_hess)
#print stepNewton(propane,xs,propane_der, propane_hess)
#print DFP(propane, x0, propane_der, propane_hess, maxiter = 2000)
#print stepNewton(propane, x0, propane_der, propane_hess)
#print adjustedNewton(propane,xs,propane_der, propane_hess)
#print BFGS(propane,xs,propane_der, propane_hess)
#print SR1(propane, x0, propane_der, propane_hess, 1e-14)
#res = scipy.optimize.minimize(propane, xs, method='BFGS',jac = propane_der,
                              #options={'disp':True})
#print (res.x)

