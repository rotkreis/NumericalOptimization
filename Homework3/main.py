import numpy as np
from penalty import BoundPenalty
from propane import propane, propane_der, propane_hess, propane_cons,propane_cons_der

xsol = np.array([0.31e-2 , 0.345e2, 0.65e-1, .859, .369e-1])
xs = 10 * xsol
#propane(xs)
#propane_der(xs)
#propane_hess(xs)
res = BoundPenalty(propane,xs,propane_der, propane_cons, propane_cons_der)
#print propane_cons_hess(xs)
print res


