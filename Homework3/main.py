import numpy as np
from penalty import BoundPenalty, Lagrange, ExPenalty
import thermo
import propane
import enzyme

"""Propane"""
xsol = np.array([0.31e-2 , 0.345e2, 0.65e-1, .859, .369e-1])
xs = 20 * xsol
#res = BoundPenalty(propane.f,xs,propane.fprime, propane.cons, propane.cons_der, wolfe2=1e-1)
#res = Lagrange(propane.f, xs, propane.fprime, propane.ins,
               #sigma0 = 1, gamma0 = np.ones(len(propane.ins(xs))),
               #diag = False, wolfe2 = .8)
#res = ExPenalty(propane.f, xs, propane.fprime, propane.ins, propane.ins_der, diag = False)

"""Thermo"""
xs = np.array([2.0e-2, 4.0e2, 2.5e2])
res = BoundPenalty(thermo.f, xs, thermo.fprime, thermo.cons, thermo.cons_der)
#res = Lagrange(thermo.f, xs, thermo.fprime, thermo.ins,
               #sigma0 = 1, gamma0 = np.ones(len(thermo.ins(xs))),
               #diag = True, wolfe2 =.91)
#res = ExPenalty(thermo.f, xs, thermo.fprime, thermo.ins, thermo.ins_der, wolfe2 = 4e-1)
#for i in range(0,len(res)):
    #print res[i]


"""Enzyme"""
#xs = np.array([2.5e-1, 3.9e-1, 4.15e-1, 3.9e-1])
#xs = 50*xs
#res = BoundPenalty(enzyme.f, xs, enzyme.fprime, enzyme.cons, enzyme.cons_der)
#res = ExPenalty(enzyme.f, xs, enzyme.fprime, enzyme.ins, enzyme.ins_der, diag = True)
#res = Lagrange(enzyme.f, xs, enzyme.fprime, enzyme.ins,
               #sigma0 = 1, gamma0 = np.ones(len(enzyme.ins(xs))),
               #diag = True)
##for i in range(0,len(res)):
    #print res[i]


