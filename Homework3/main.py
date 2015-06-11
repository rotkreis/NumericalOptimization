import numpy as np
from penalty import BoundPenalty, Lagrange, ExPenalty
from propane import propane, propane_der, propane_hess, propane_cons,propane_cons_der
import thermo
import enzyme

"""Propane"""
#xsol = np.array([0.31e-2 , 0.345e2, 0.65e-1, .859, .369e-1])
#xs = 21 * xsol
#propane(xs)
#propane_der(xs)
#propane_hess(xs)
#res = BoundPenalty(propane,xs,propane_der, propane_cons, propane_cons_der)
#res = Lagrange(propane, xs, propane_der, propane_cons,
               #sigma0 = 1, gamma0 = np.ones(len(propane_cons(xs))),
               #diag = True)

"""Thermo"""
#xs = np.array([2.0e-2, 4.0e2, 2.5e2])
##res = BoundPenalty(thermo.f, xs, thermo.fprime, thermo.cons, thermo.cons_der)
#res = Lagrange(thermo.f, xs, thermo.fprime, thermo.ins,
               #sigma0 = 1, gamma0 = np.ones(len(thermo.ins(xs))),
               #diag = True)

#for i in range(0,len(res)):
    #print res[i]


"""Enzyme"""
xs = np.array([2.5e-1, 3.9e-1, 4.15e-1, 3.9e-1])
xs = 50*xs
#res = BoundPenalty(enzyme.f, xs, enzyme.fprime, enzyme.cons, enzyme.cons_der)
res = ExPenalty(enzyme.f, xs, enzyme.fprime, enzyme.ins, enzyme.ins_der)
#res = Lagrange(enzyme.f, xs, enzyme.fprime, enzyme.ins,
               #sigma0 = 1, gamma0 = np.ones(len(enzyme.ins(xs))),
               #diag = True)
#for i in range(0,len(res)):
    #print res[i]


