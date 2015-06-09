import numpy as np
from penalty import BoundPenalty, Lagrange, ExPenalty
from propane import propane, propane_der, propane_hess, propane_cons,propane_cons_der

xsol = np.array([0.31e-2 , 0.345e2, 0.65e-1, .859, .369e-1])
xs = 21 * xsol
#propane(xs)
#propane_der(xs)
#propane_hess(xs)
#res = BoundPenalty(propane,xs,propane_der, propane_cons, propane_cons_der)
res = Lagrange(propane, xs, propane_der, propane_cons,
               sigma0 = 1, gamma0 = np.ones(len(propane_cons(xs))),
               diag = True)
#print propane_cons_hess(xs)
for i in range(0,len(res)):
    print res[i]


