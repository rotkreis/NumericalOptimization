import numpy as np
sigma1=np.array([.485, -.0019, -.0581, .015, .105, .0406, .167, -.399])
x1=np.array([.299, .186, -.0273, .0254, -.474, .474, -.892, 0.892])
sigma2=np.array([-.69, -.044, -1.57, -1.31, -2.65, 2.0, -12.6, 9.48])
x2=np.array([-.3, -.39, .3, -3.44, -1.2, 2.69, 1.59, -1.5])
sigma3=np.array([-.816, -.017, -1.826, -.754,-4.839, -3.259,-14.023,15.467])
x3=np.array([-.041,-.775,.03, -.047, -2.565, 2.565, -.754, .754])
sig_x,sig_y, sig_a,sig_b,sig_c,sig_d,sig_e,sig_f = sigma3
def dipole(x):
    x1 = x[0]
    x2 = x[1]
    x3 = x[2]
    x4 = x[3]
    x5 = x[4]
    x6 = x[5]
    x7 = x[6]
    x8 = x[7]
    f1 = x1 + x2 - sig_x
    f2 = x3 + x4 - sig_y
    f3 = x5*x1 + x6*x2 - x7*x3 - x8*x4 - sig_a
    f4 = x7*x1 + x8*x2 + x5*x3 + x6*x4 - sig_b
    f5 = x1*(x5**2-x7**2) - 2*x3*x5*x7 + x2*(x6**2-x8**2) - 2*x4*x6*x8 - sig_c
    f6 = x3*(x5**2-x7**2) + 2*x1*x5*x7 + x4*(x6**2-x8**2) + 2*x2*x6*x8 - sig_d
    f7 = x1*x5*(x5**2 - 3*x7**2) + x3*x7*(x7**2 - 3*x5**2) + x2*x6*(x6**2 - 3*x8**2) + x4*x8*(x8**2 - 3*x6**2) - sig_e
    f8 = x3*x5*(x5**2 - 3*x7**2) - x1*x7*(x7**2 - 3*x5**2) + x4*x6*(x6**2 - 3*x8**2) - x2*x8*(x8**2 - 3*x6**2) - sig_f
    return f1**2 + f2**2 + f3**2 + f4**2 + f5**2 + f6**2 + f7**2 + f8**2
def dipole_der(x):
    x1 = x[0]
    x2 = x[1]
    x3 = x[2]
    x4 = x[3]
    x5 = x[4]
    x6 = x[5]
    x7 = x[6]
    x8 = x[7]
    der = np.zeros_like(x)
    der[0] =  -2*sig_x + 2*x1 + 2*x2 + 4*x5*x7*(-sig_d + 2*x1*x5*x7 + 2*x2*x6*x8 + x3*(x5**2 - x7**2) + x4*(x6**2 - x8**2)) + 2*x5*(x5**2 - 3*x7**2)*(-sig_e + x1*x5*(x5**2 - 3*x7**2) + x2*x6*(x6**2 - 3*x8**2) + x3*x7*(-3*x5**2 + x7**2) + x4*x8*(-3*x6**2 + x8**2)) + 2*x5*(-sig_a + x1*x5 + x2*x6 - x3*x7 - x4*x8) - 2*x7*(-3*x5**2 + x7**2)*(-sig_f - x1*x7*(-3*x5**2 + x7**2) - x2*x8*(-3*x6**2 + x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2)) + 2*x7*(-sig_b + x1*x7 + x2*x8 + x3*x5 + x4*x6) + (2*x5**2 - 2*x7**2)*(-sig_c + x1*(x5**2 - x7**2) + x2*(x6**2 - x8**2) - 2*x3*x5*x7 - 2*x4*x6*x8)
    der[1] =  -2*sig_x + 2*x1 + 2*x2 + 4*x6*x8*(-sig_d + 2*x1*x5*x7 + 2*x2*x6*x8 + x3*(x5**2 - x7**2) + x4*(x6**2 - x8**2)) + 2*x6*(x6**2 - 3*x8**2)*(-sig_e + x1*x5*(x5**2 - 3*x7**2) + x2*x6*(x6**2 - 3*x8**2) + x3*x7*(-3*x5**2 + x7**2) + x4*x8*(-3*x6**2 + x8**2)) + 2*x6*(-sig_a + x1*x5 + x2*x6 - x3*x7 - x4*x8) - 2*x8*(-3*x6**2 + x8**2)*(-sig_f - x1*x7*(-3*x5**2 + x7**2) - x2*x8*(-3*x6**2 + x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2)) + 2*x8*(-sig_b + x1*x7 + x2*x8 + x3*x5 + x4*x6) + (2*x6**2 - 2*x8**2)*(-sig_c + x1*(x5**2 - x7**2) + x2*(x6**2 - x8**2) - 2*x3*x5*x7 - 2*x4*x6*x8)
    der[2] =  -2*sig_y + 2*x3 + 2*x4 - 4*x5*x7*(-sig_c + x1*(x5**2 - x7**2) + x2*(x6**2 - x8**2) - 2*x3*x5*x7 - 2*x4*x6*x8) + 2*x5*(x5**2 - 3*x7**2)*(-sig_f - x1*x7*(-3*x5**2 + x7**2) - x2*x8*(-3*x6**2 + x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2)) + 2*x5*(-sig_b + x1*x7 + x2*x8 + x3*x5 + x4*x6) + 2*x7*(-3*x5**2 + x7**2)*(-sig_e + x1*x5*(x5**2 - 3*x7**2) + x2*x6*(x6**2 - 3*x8**2) + x3*x7*(-3*x5**2 + x7**2) + x4*x8*(-3*x6**2 + x8**2)) - 2*x7*(-sig_a + x1*x5 + x2*x6 - x3*x7 - x4*x8) + (2*x5**2 - 2*x7**2)*(-sig_d + 2*x1*x5*x7 + 2*x2*x6*x8 + x3*(x5**2 - x7**2) + x4*(x6**2 - x8**2))
    der[3] =  -2*sig_y + 2*x3 + 2*x4 - 4*x6*x8*(-sig_c + x1*(x5**2 - x7**2) + x2*(x6**2 - x8**2) - 2*x3*x5*x7 - 2*x4*x6*x8) + 2*x6*(x6**2 - 3*x8**2)*(-sig_f - x1*x7*(-3*x5**2 + x7**2) - x2*x8*(-3*x6**2 + x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2)) + 2*x6*(-sig_b + x1*x7 + x2*x8 + x3*x5 + x4*x6) + 2*x8*(-3*x6**2 + x8**2)*(-sig_e + x1*x5*(x5**2 - 3*x7**2) + x2*x6*(x6**2 - 3*x8**2) + x3*x7*(-3*x5**2 + x7**2) + x4*x8*(-3*x6**2 + x8**2)) - 2*x8*(-sig_a + x1*x5 + x2*x6 - x3*x7 - x4*x8) + (2*x6**2 - 2*x8**2)*(-sig_d + 2*x1*x5*x7 + 2*x2*x6*x8 + x3*(x5**2 - x7**2) + x4*(x6**2 - x8**2))
    der[4] =  2*x1*(-sig_a + x1*x5 + x2*x6 - x3*x7 - x4*x8) + 2*x3*(-sig_b + x1*x7 + x2*x8 + x3*x5 + x4*x6) + (4*x1*x5 - 4*x3*x7)*(-sig_c + x1*(x5**2 - x7**2) + x2*(x6**2 - x8**2) - 2*x3*x5*x7 - 2*x4*x6*x8) + (4*x1*x7 + 4*x3*x5)*(-sig_d + 2*x1*x5*x7 + 2*x2*x6*x8 + x3*(x5**2 - x7**2) + x4*(x6**2 - x8**2)) + (4*x1*x5**2 + 2*x1*(x5**2 - 3*x7**2) - 12*x3*x5*x7)*(-sig_e + x1*x5*(x5**2 - 3*x7**2) + x2*x6*(x6**2 - 3*x8**2) + x3*x7*(-3*x5**2 + x7**2) + x4*x8*(-3*x6**2 + x8**2)) + (12*x1*x5*x7 + 4*x3*x5**2 + 2*x3*(x5**2 - 3*x7**2))*(-sig_f - x1*x7*(-3*x5**2 + x7**2) - x2*x8*(-3*x6**2 + x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2))
    der[5] =  2*x2*(-sig_a + x1*x5 + x2*x6 - x3*x7 - x4*x8) + 2*x4*(-sig_b + x1*x7 + x2*x8 + x3*x5 + x4*x6) + (4*x2*x6 - 4*x4*x8)*(-sig_c + x1*(x5**2 - x7**2) + x2*(x6**2 - x8**2) - 2*x3*x5*x7 - 2*x4*x6*x8) + (4*x2*x8 + 4*x4*x6)*(-sig_d + 2*x1*x5*x7 + 2*x2*x6*x8 + x3*(x5**2 - x7**2) + x4*(x6**2 - x8**2)) + (4*x2*x6**2 + 2*x2*(x6**2 - 3*x8**2) - 12*x4*x6*x8)*(-sig_e + x1*x5*(x5**2 - 3*x7**2) + x2*x6*(x6**2 - 3*x8**2) + x3*x7*(-3*x5**2 + x7**2) + x4*x8*(-3*x6**2 + x8**2)) + (12*x2*x6*x8 + 4*x4*x6**2 + 2*x4*(x6**2 - 3*x8**2))*(-sig_f - x1*x7*(-3*x5**2 + x7**2) - x2*x8*(-3*x6**2 + x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2))
    der[6] =  2*x1*(-sig_b + x1*x7 + x2*x8 + x3*x5 + x4*x6) - 2*x3*(-sig_a + x1*x5 + x2*x6 - x3*x7 - x4*x8) + (4*x1*x5 - 4*x3*x7)*(-sig_d + 2*x1*x5*x7 + 2*x2*x6*x8 + x3*(x5**2 - x7**2) + x4*(x6**2 - x8**2)) + (-4*x1*x7 - 4*x3*x5)*(-sig_c + x1*(x5**2 - x7**2) + x2*(x6**2 - x8**2) - 2*x3*x5*x7 - 2*x4*x6*x8) + (-4*x1*x7**2 - 2*x1*(-3*x5**2 + x7**2) - 12*x3*x5*x7)*(-sig_f - x1*x7*(-3*x5**2 + x7**2) - x2*x8*(-3*x6**2 + x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2)) + (-12*x1*x5*x7 + 4*x3*x7**2 + 2*x3*(-3*x5**2 + x7**2))*(-sig_e + x1*x5*(x5**2 - 3*x7**2) + x2*x6*(x6**2 - 3*x8**2) + x3*x7*(-3*x5**2 + x7**2) + x4*x8*(-3*x6**2 + x8**2))
    der[7] =  2*x2*(-sig_b + x1*x7 + x2*x8 + x3*x5 + x4*x6) - 2*x4*(-sig_a + x1*x5 + x2*x6 - x3*x7 - x4*x8) + (4*x2*x6 - 4*x4*x8)*(-sig_d + 2*x1*x5*x7 + 2*x2*x6*x8 + x3*(x5**2 - x7**2) + x4*(x6**2 - x8**2)) + (-4*x2*x8 - 4*x4*x6)*(-sig_c + x1*(x5**2 - x7**2) + x2*(x6**2 - x8**2) - 2*x3*x5*x7 - 2*x4*x6*x8) + (-4*x2*x8**2 - 2*x2*(-3*x6**2 + x8**2) - 12*x4*x6*x8)*(-sig_f - x1*x7*(-3*x5**2 + x7**2) - x2*x8*(-3*x6**2 + x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2)) + (-12*x2*x6*x8 + 4*x4*x8**2 + 2*x4*(-3*x6**2 + x8**2))*(-sig_e + x1*x5*(x5**2 - 3*x7**2) + x2*x6*(x6**2 - 3*x8**2) + x3*x7*(-3*x5**2 + x7**2) + x4*x8*(-3*x6**2 + x8**2))
    return der
def dipole_hess(x):
    x1 = x[0]
    x2 = x[1]
    x3 = x[2]
    x4 = x[3]
    x5 = x[4]
    x6 = x[5]
    x7 = x[6]
    x8 = x[7]
    hess = np.zeros((8,8))
    hess[0][0] =  2*(4*x5**2*x7**2 + x5**2*(x5**2 - 3*x7**2)**2 + x5**2 + x7**2*(3*x5**2 - x7**2)**2 + x7**2 + (x5**2 - x7**2)**2 + 1)
    hess[0][1] =  2*(4*x5*x6*x7*x8 + x5*x6*(x5**2 - 3*x7**2)*(x6**2 - 3*x8**2) + x5*x6 + x7*x8*(3*x5**2 - x7**2)*(3*x6**2 - x8**2) + x7*x8 + (x5**2 - x7**2)*(x6**2 - x8**2) + 1)
    hess[0][2] =  0
    hess[0][3] =  2*(2*x5*x7*(x6**2 - x8**2) - x5*x8*(x5**2 - 3*x7**2)*(3*x6**2 - x8**2) - x5*x8 + x6*x7*(3*x5**2 - x7**2)*(x6**2 - 3*x8**2) + x6*x7 - 2*x6*x8*(x5**2 - x7**2))
    hess[0][4] =  2*(-sig_a + 2*x1*x5 + x2*x6 - x4*x8 - 2*x5**2*(sig_e - x1*x5*(x5**2 - 3*x7**2) - x2*x6*(x6**2 - 3*x8**2) + x3*x7*(3*x5**2 - x7**2) + x4*x8*(3*x6**2 - x8**2)) + 4*x5*x7*(x1*x7 + x3*x5) + 6*x5*x7*(-sig_f + x1*x7*(3*x5**2 - x7**2) + x2*x8*(3*x6**2 - x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2)) + x5*(x5**2 - 3*x7**2)*(2*x1*x5**2 + x1*(x5**2 - 3*x7**2) - 6*x3*x5*x7) - 2*x5*(sig_c - x1*(x5**2 - x7**2) - x2*(x6**2 - x8**2) + 2*x3*x5*x7 + 2*x4*x6*x8) + x7*(3*x5**2 - x7**2)*(6*x1*x5*x7 + 2*x3*x5**2 + x3*(x5**2 - 3*x7**2)) + 2*x7*(-sig_d + 2*x1*x5*x7 + 2*x2*x6*x8 + x3*(x5**2 - x7**2) + x4*(x6**2 - x8**2)) - (x5**2 - 3*x7**2)*(sig_e - x1*x5*(x5**2 - 3*x7**2) - x2*x6*(x6**2 - 3*x8**2) + x3*x7*(3*x5**2 - x7**2) + x4*x8*(3*x6**2 - x8**2)) + 2*(x5**2 - x7**2)*(x1*x5 - x3*x7))
    hess[0][5] =  2*(x2*x5 + x4*x7 + 4*x5*x7*(x2*x8 + x4*x6) + x5*(x5**2 - 3*x7**2)*(2*x2*x6**2 + x2*(x6**2 - 3*x8**2) - 6*x4*x6*x8) + x7*(3*x5**2 - x7**2)*(6*x2*x6*x8 + 2*x4*x6**2 + x4*(x6**2 - 3*x8**2)) + 2*(x5**2 - x7**2)*(x2*x6 - x4*x8))
    hess[0][6] =  2*(-sig_b + 2*x1*x7 + x2*x8 + x4*x6 + 4*x5*x7*(x1*x5 - x3*x7) + 6*x5*x7*(sig_e - x1*x5*(x5**2 - 3*x7**2) - x2*x6*(x6**2 - 3*x8**2) + x3*x7*(3*x5**2 - x7**2) + x4*x8*(3*x6**2 - x8**2)) - x5*(x5**2 - 3*x7**2)*(6*x1*x5*x7 - 2*x3*x7**2 + x3*(3*x5**2 - x7**2)) + 2*x5*(-sig_d + 2*x1*x5*x7 + 2*x2*x6*x8 + x3*(x5**2 - x7**2) + x4*(x6**2 - x8**2)) - 2*x7**2*(-sig_f + x1*x7*(3*x5**2 - x7**2) + x2*x8*(3*x6**2 - x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2)) - x7*(3*x5**2 - x7**2)*(2*x1*x7**2 - x1*(3*x5**2 - x7**2) + 6*x3*x5*x7) + 2*x7*(sig_c - x1*(x5**2 - x7**2) - x2*(x6**2 - x8**2) + 2*x3*x5*x7 + 2*x4*x6*x8) - 2*(x5**2 - x7**2)*(x1*x7 + x3*x5) + (3*x5**2 - x7**2)*(-sig_f + x1*x7*(3*x5**2 - x7**2) + x2*x8*(3*x6**2 - x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2)))
    hess[0][7] =  2*(x2*x7 - x4*x5 + 4*x5*x7*(x2*x6 - x4*x8) - x5*(x5**2 - 3*x7**2)*(6*x2*x6*x8 - 2*x4*x8**2 + x4*(3*x6**2 - x8**2)) - x7*(3*x5**2 - x7**2)*(2*x2*x8**2 - x2*(3*x6**2 - x8**2) + 6*x4*x6*x8) - 2*(x5**2 - x7**2)*(x2*x8 + x4*x6))
    hess[1][0] =  2*(4*x5*x6*x7*x8 + x5*x6*(x5**2 - 3*x7**2)*(x6**2 - 3*x8**2) + x5*x6 + x7*x8*(3*x5**2 - x7**2)*(3*x6**2 - x8**2) + x7*x8 + (x5**2 - x7**2)*(x6**2 - x8**2) + 1)
    hess[1][1] =  2*(4*x6**2*x8**2 + x6**2*(x6**2 - 3*x8**2)**2 + x6**2 + x8**2*(3*x6**2 - x8**2)**2 + x8**2 + (x6**2 - x8**2)**2 + 1)
    hess[1][2] =  2*(-2*x5*x7*(x6**2 - x8**2) + x5*x8*(x5**2 - 3*x7**2)*(3*x6**2 - x8**2) + x5*x8 - x6*x7*(3*x5**2 - x7**2)*(x6**2 - 3*x8**2) - x6*x7 + 2*x6*x8*(x5**2 - x7**2))
    hess[1][3] =  0
    hess[1][4] =  2*(x1*x6 + x3*x8 + 4*x6*x8*(x1*x7 + x3*x5) + x6*(x6**2 - 3*x8**2)*(2*x1*x5**2 + x1*(x5**2 - 3*x7**2) - 6*x3*x5*x7) + x8*(3*x6**2 - x8**2)*(6*x1*x5*x7 + 2*x3*x5**2 + x3*(x5**2 - 3*x7**2)) + 2*(x6**2 - x8**2)*(x1*x5 - x3*x7))
    hess[1][5] =  2*(-sig_a + x1*x5 + 2*x2*x6 - x3*x7 - 2*x6**2*(sig_e - x1*x5*(x5**2 - 3*x7**2) - x2*x6*(x6**2 - 3*x8**2) + x3*x7*(3*x5**2 - x7**2) + x4*x8*(3*x6**2 - x8**2)) + 4*x6*x8*(x2*x8 + x4*x6) + 6*x6*x8*(-sig_f + x1*x7*(3*x5**2 - x7**2) + x2*x8*(3*x6**2 - x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2)) + x6*(x6**2 - 3*x8**2)*(2*x2*x6**2 + x2*(x6**2 - 3*x8**2) - 6*x4*x6*x8) - 2*x6*(sig_c - x1*(x5**2 - x7**2) - x2*(x6**2 - x8**2) + 2*x3*x5*x7 + 2*x4*x6*x8) + x8*(3*x6**2 - x8**2)*(6*x2*x6*x8 + 2*x4*x6**2 + x4*(x6**2 - 3*x8**2)) + 2*x8*(-sig_d + 2*x1*x5*x7 + 2*x2*x6*x8 + x3*(x5**2 - x7**2) + x4*(x6**2 - x8**2)) - (x6**2 - 3*x8**2)*(sig_e - x1*x5*(x5**2 - 3*x7**2) - x2*x6*(x6**2 - 3*x8**2) + x3*x7*(3*x5**2 - x7**2) + x4*x8*(3*x6**2 - x8**2)) + 2*(x6**2 - x8**2)*(x2*x6 - x4*x8))
    hess[1][6] =  2*(x1*x8 - x3*x6 + 4*x6*x8*(x1*x5 - x3*x7) - x6*(x6**2 - 3*x8**2)*(6*x1*x5*x7 - 2*x3*x7**2 + x3*(3*x5**2 - x7**2)) - x8*(3*x6**2 - x8**2)*(2*x1*x7**2 - x1*(3*x5**2 - x7**2) + 6*x3*x5*x7) - 2*(x6**2 - x8**2)*(x1*x7 + x3*x5))
    hess[1][7] =  2*(-sig_b + x1*x7 + 2*x2*x8 + x3*x5 + 4*x6*x8*(x2*x6 - x4*x8) + 6*x6*x8*(sig_e - x1*x5*(x5**2 - 3*x7**2) - x2*x6*(x6**2 - 3*x8**2) + x3*x7*(3*x5**2 - x7**2) + x4*x8*(3*x6**2 - x8**2)) - x6*(x6**2 - 3*x8**2)*(6*x2*x6*x8 - 2*x4*x8**2 + x4*(3*x6**2 - x8**2)) + 2*x6*(-sig_d + 2*x1*x5*x7 + 2*x2*x6*x8 + x3*(x5**2 - x7**2) + x4*(x6**2 - x8**2)) - 2*x8**2*(-sig_f + x1*x7*(3*x5**2 - x7**2) + x2*x8*(3*x6**2 - x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2)) - x8*(3*x6**2 - x8**2)*(2*x2*x8**2 - x2*(3*x6**2 - x8**2) + 6*x4*x6*x8) + 2*x8*(sig_c - x1*(x5**2 - x7**2) - x2*(x6**2 - x8**2) + 2*x3*x5*x7 + 2*x4*x6*x8) - 2*(x6**2 - x8**2)*(x2*x8 + x4*x6) + (3*x6**2 - x8**2)*(-sig_f + x1*x7*(3*x5**2 - x7**2) + x2*x8*(3*x6**2 - x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2)))
    hess[2][0] =  0
    hess[2][1] =  2*(-2*x5*x7*(x6**2 - x8**2) + x5*x8*(x5**2 - 3*x7**2)*(3*x6**2 - x8**2) + x5*x8 - x6*x7*(3*x5**2 - x7**2)*(x6**2 - 3*x8**2) - x6*x7 + 2*x6*x8*(x5**2 - x7**2))
    hess[2][2] =  2*(4*x5**2*x7**2 + x5**2*(x5**2 - 3*x7**2)**2 + x5**2 + x7**2*(3*x5**2 - x7**2)**2 + x7**2 + (x5**2 - x7**2)**2 + 1)
    hess[2][3] =  2*(4*x5*x6*x7*x8 + x5*x6*(x5**2 - 3*x7**2)*(x6**2 - 3*x8**2) + x5*x6 + x7*x8*(3*x5**2 - x7**2)*(3*x6**2 - x8**2) + x7*x8 + (x5**2 - x7**2)*(x6**2 - x8**2) + 1)
    hess[2][4] =  2*(-sig_b + x2*x8 + 2*x3*x5 + x4*x6 + 2*x5**2*(-sig_f + x1*x7*(3*x5**2 - x7**2) + x2*x8*(3*x6**2 - x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2)) - 4*x5*x7*(x1*x5 - x3*x7) + 6*x5*x7*(sig_e - x1*x5*(x5**2 - 3*x7**2) - x2*x6*(x6**2 - 3*x8**2) + x3*x7*(3*x5**2 - x7**2) + x4*x8*(3*x6**2 - x8**2)) + x5*(x5**2 - 3*x7**2)*(6*x1*x5*x7 + 2*x3*x5**2 + x3*(x5**2 - 3*x7**2)) + 2*x5*(-sig_d + 2*x1*x5*x7 + 2*x2*x6*x8 + x3*(x5**2 - x7**2) + x4*(x6**2 - x8**2)) - x7*(3*x5**2 - x7**2)*(2*x1*x5**2 + x1*(x5**2 - 3*x7**2) - 6*x3*x5*x7) + 2*x7*(sig_c - x1*(x5**2 - x7**2) - x2*(x6**2 - x8**2) + 2*x3*x5*x7 + 2*x4*x6*x8) + (x5**2 - 3*x7**2)*(-sig_f + x1*x7*(3*x5**2 - x7**2) + x2*x8*(3*x6**2 - x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2)) + 2*(x5**2 - x7**2)*(x1*x7 + x3*x5))
    hess[2][5] =  2*(-x2*x7 + x4*x5 - 4*x5*x7*(x2*x6 - x4*x8) + x5*(x5**2 - 3*x7**2)*(6*x2*x6*x8 + 2*x4*x6**2 + x4*(x6**2 - 3*x8**2)) - x7*(3*x5**2 - x7**2)*(2*x2*x6**2 + x2*(x6**2 - 3*x8**2) - 6*x4*x6*x8) + 2*(x5**2 - x7**2)*(x2*x8 + x4*x6))
    hess[2][6] =  2*(sig_a - x2*x6 + 2*x3*x7 + x4*x8 + 4*x5*x7*(x1*x7 + x3*x5) - 6*x5*x7*(-sig_f + x1*x7*(3*x5**2 - x7**2) + x2*x8*(3*x6**2 - x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2)) - x5*(x5**2 - 3*x7**2)*(2*x1*x7**2 - x1*(3*x5**2 - x7**2) + 6*x3*x5*x7) + 2*x5*(sig_c - x1*(x5**2 - x7**2) - x2*(x6**2 - x8**2) + 2*x3*x5*x7 + 2*x4*x6*x8) - 2*x7**2*(sig_e - x1*x5*(x5**2 - 3*x7**2) - x2*x6*(x6**2 - 3*x8**2) + x3*x7*(3*x5**2 - x7**2) + x4*x8*(3*x6**2 - x8**2)) + x7*(3*x5**2 - x7**2)*(6*x1*x5*x7 - 2*x3*x7**2 + x3*(3*x5**2 - x7**2)) - 2*x7*(-sig_d + 2*x1*x5*x7 + 2*x2*x6*x8 + x3*(x5**2 - x7**2) + x4*(x6**2 - x8**2)) + 2*(x5**2 - x7**2)*(x1*x5 - x3*x7) + (3*x5**2 - x7**2)*(sig_e - x1*x5*(x5**2 - 3*x7**2) - x2*x6*(x6**2 - 3*x8**2) + x3*x7*(3*x5**2 - x7**2) + x4*x8*(3*x6**2 - x8**2)))
    hess[2][7] =  2*(x2*x5 + x4*x7 + 4*x5*x7*(x2*x8 + x4*x6) - x5*(x5**2 - 3*x7**2)*(2*x2*x8**2 - x2*(3*x6**2 - x8**2) + 6*x4*x6*x8) + x7*(3*x5**2 - x7**2)*(6*x2*x6*x8 - 2*x4*x8**2 + x4*(3*x6**2 - x8**2)) + 2*(x5**2 - x7**2)*(x2*x6 - x4*x8))
    hess[3][0] =  2*(2*x5*x7*(x6**2 - x8**2) - x5*x8*(x5**2 - 3*x7**2)*(3*x6**2 - x8**2) - x5*x8 + x6*x7*(3*x5**2 - x7**2)*(x6**2 - 3*x8**2) + x6*x7 - 2*x6*x8*(x5**2 - x7**2))
    hess[3][1] =  0
    hess[3][2] =  2*(4*x5*x6*x7*x8 + x5*x6*(x5**2 - 3*x7**2)*(x6**2 - 3*x8**2) + x5*x6 + x7*x8*(3*x5**2 - x7**2)*(3*x6**2 - x8**2) + x7*x8 + (x5**2 - x7**2)*(x6**2 - x8**2) + 1)
    hess[3][3] =  2*(4*x6**2*x8**2 + x6**2*(x6**2 - 3*x8**2)**2 + x6**2 + x8**2*(3*x6**2 - x8**2)**2 + x8**2 + (x6**2 - x8**2)**2 + 1)
    hess[3][4] =  2*(-x1*x8 + x3*x6 - 4*x6*x8*(x1*x5 - x3*x7) + x6*(x6**2 - 3*x8**2)*(6*x1*x5*x7 + 2*x3*x5**2 + x3*(x5**2 - 3*x7**2)) - x8*(3*x6**2 - x8**2)*(2*x1*x5**2 + x1*(x5**2 - 3*x7**2) - 6*x3*x5*x7) + 2*(x6**2 - x8**2)*(x1*x7 + x3*x5))
    hess[3][5] =  2*(-sig_b + x1*x7 + x3*x5 + 2*x4*x6 + 2*x6**2*(-sig_f + x1*x7*(3*x5**2 - x7**2) + x2*x8*(3*x6**2 - x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2)) - 4*x6*x8*(x2*x6 - x4*x8) + 6*x6*x8*(sig_e - x1*x5*(x5**2 - 3*x7**2) - x2*x6*(x6**2 - 3*x8**2) + x3*x7*(3*x5**2 - x7**2) + x4*x8*(3*x6**2 - x8**2)) + x6*(x6**2 - 3*x8**2)*(6*x2*x6*x8 + 2*x4*x6**2 + x4*(x6**2 - 3*x8**2)) + 2*x6*(-sig_d + 2*x1*x5*x7 + 2*x2*x6*x8 + x3*(x5**2 - x7**2) + x4*(x6**2 - x8**2)) - x8*(3*x6**2 - x8**2)*(2*x2*x6**2 + x2*(x6**2 - 3*x8**2) - 6*x4*x6*x8) + 2*x8*(sig_c - x1*(x5**2 - x7**2) - x2*(x6**2 - x8**2) + 2*x3*x5*x7 + 2*x4*x6*x8) + (x6**2 - 3*x8**2)*(-sig_f + x1*x7*(3*x5**2 - x7**2) + x2*x8*(3*x6**2 - x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2)) + 2*(x6**2 - x8**2)*(x2*x8 + x4*x6))
    hess[3][6] =  2*(x1*x6 + x3*x8 + 4*x6*x8*(x1*x7 + x3*x5) - x6*(x6**2 - 3*x8**2)*(2*x1*x7**2 - x1*(3*x5**2 - x7**2) + 6*x3*x5*x7) + x8*(3*x6**2 - x8**2)*(6*x1*x5*x7 - 2*x3*x7**2 + x3*(3*x5**2 - x7**2)) + 2*(x6**2 - x8**2)*(x1*x5 - x3*x7))
    hess[3][7] =  2*(sig_a - x1*x5 + x3*x7 + 2*x4*x8 + 4*x6*x8*(x2*x8 + x4*x6) - 6*x6*x8*(-sig_f + x1*x7*(3*x5**2 - x7**2) + x2*x8*(3*x6**2 - x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2)) - x6*(x6**2 - 3*x8**2)*(2*x2*x8**2 - x2*(3*x6**2 - x8**2) + 6*x4*x6*x8) + 2*x6*(sig_c - x1*(x5**2 - x7**2) - x2*(x6**2 - x8**2) + 2*x3*x5*x7 + 2*x4*x6*x8) - 2*x8**2*(sig_e - x1*x5*(x5**2 - 3*x7**2) - x2*x6*(x6**2 - 3*x8**2) + x3*x7*(3*x5**2 - x7**2) + x4*x8*(3*x6**2 - x8**2)) + x8*(3*x6**2 - x8**2)*(6*x2*x6*x8 - 2*x4*x8**2 + x4*(3*x6**2 - x8**2)) - 2*x8*(-sig_d + 2*x1*x5*x7 + 2*x2*x6*x8 + x3*(x5**2 - x7**2) + x4*(x6**2 - x8**2)) + 2*(x6**2 - x8**2)*(x2*x6 - x4*x8) + (3*x6**2 - x8**2)*(sig_e - x1*x5*(x5**2 - 3*x7**2) - x2*x6*(x6**2 - 3*x8**2) + x3*x7*(3*x5**2 - x7**2) + x4*x8*(3*x6**2 - x8**2)))
    hess[4][0] =  2*(-sig_a + 2*x1*x5 + x2*x6 - x4*x8 + 4*x5*x7*(x1*x7 + x3*x5) + 6*x5*x7*(-sig_f + x1*x7*(3*x5**2 - x7**2) + x2*x8*(3*x6**2 - x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2)) + x5*(x5**2 - 3*x7**2)*(2*x1*x5**2 + x1*(x5**2 - 3*x7**2) - 6*x3*x5*x7) - 2*x5*(sig_c - x1*(x5**2 - x7**2) - x2*(x6**2 - x8**2) + 2*x3*x5*x7 + 2*x4*x6*x8) + x7*(3*x5**2 - x7**2)*(6*x1*x5*x7 + 2*x3*x5**2 + x3*(x5**2 - 3*x7**2)) + 2*x7*(-sig_d + 2*x1*x5*x7 + 2*x2*x6*x8 + x3*(x5**2 - x7**2) + x4*(x6**2 - x8**2)) + 2*(x5**2 - x7**2)*(x1*x5 - x3*x7) - 3*(x5**2 - x7**2)*(sig_e - x1*x5*(x5**2 - 3*x7**2) - x2*x6*(x6**2 - 3*x8**2) + x3*x7*(3*x5**2 - x7**2) + x4*x8*(3*x6**2 - x8**2)))
    hess[4][1] =  2*(x1*x6 + x3*x8 + 4*x6*x8*(x1*x7 + x3*x5) + x6*(x6**2 - 3*x8**2)*(2*x1*x5**2 + x1*(x5**2 - 3*x7**2) - 6*x3*x5*x7) + x8*(3*x6**2 - x8**2)*(6*x1*x5*x7 + 2*x3*x5**2 + x3*(x5**2 - 3*x7**2)) + 2*(x6**2 - x8**2)*(x1*x5 - x3*x7))
    hess[4][2] =  2*(-sig_b + x2*x8 + 2*x3*x5 + x4*x6 - 4*x5*x7*(x1*x5 - x3*x7) + 6*x5*x7*(sig_e - x1*x5*(x5**2 - 3*x7**2) - x2*x6*(x6**2 - 3*x8**2) + x3*x7*(3*x5**2 - x7**2) + x4*x8*(3*x6**2 - x8**2)) + x5*(x5**2 - 3*x7**2)*(6*x1*x5*x7 + 2*x3*x5**2 + x3*(x5**2 - 3*x7**2)) + 2*x5*(-sig_d + 2*x1*x5*x7 + 2*x2*x6*x8 + x3*(x5**2 - x7**2) + x4*(x6**2 - x8**2)) - x7*(3*x5**2 - x7**2)*(2*x1*x5**2 + x1*(x5**2 - 3*x7**2) - 6*x3*x5*x7) + 2*x7*(sig_c - x1*(x5**2 - x7**2) - x2*(x6**2 - x8**2) + 2*x3*x5*x7 + 2*x4*x6*x8) + 2*(x5**2 - x7**2)*(x1*x7 + x3*x5) + 3*(x5**2 - x7**2)*(-sig_f + x1*x7*(3*x5**2 - x7**2) + x2*x8*(3*x6**2 - x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2)))
    hess[4][3] =  2*(-x1*x8 + x3*x6 - 4*x6*x8*(x1*x5 - x3*x7) + x6*(x6**2 - 3*x8**2)*(6*x1*x5*x7 + 2*x3*x5**2 + x3*(x5**2 - 3*x7**2)) - x8*(3*x6**2 - x8**2)*(2*x1*x5**2 + x1*(x5**2 - 3*x7**2) - 6*x3*x5*x7) + 2*(x6**2 - x8**2)*(x1*x7 + x3*x5))
    hess[4][4] =  2*(x1**2 - 2*x1*(sig_c - x1*(x5**2 - x7**2) - x2*(x6**2 - x8**2) + 2*x3*x5*x7 + 2*x4*x6*x8) + x3**2 + 2*x3*(-sig_d + 2*x1*x5*x7 + 2*x2*x6*x8 + x3*(x5**2 - x7**2) + x4*(x6**2 - x8**2)) + 4*(x1*x5 - x3*x7)**2 - 6*(x1*x5 - x3*x7)*(sig_e - x1*x5*(x5**2 - 3*x7**2) - x2*x6*(x6**2 - 3*x8**2) + x3*x7*(3*x5**2 - x7**2) + x4*x8*(3*x6**2 - x8**2)) + 4*(x1*x7 + x3*x5)**2 + 6*(x1*x7 + x3*x5)*(-sig_f + x1*x7*(3*x5**2 - x7**2) + x2*x8*(3*x6**2 - x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2)) + (2*x1*x5**2 + x1*(x5**2 - 3*x7**2) - 6*x3*x5*x7)**2 + (6*x1*x5*x7 + 2*x3*x5**2 + x3*(x5**2 - 3*x7**2))**2)
    hess[4][5] =  2*(x1*x2 + x3*x4 + 4*(x1*x5 - x3*x7)*(x2*x6 - x4*x8) + 4*(x1*x7 + x3*x5)*(x2*x8 + x4*x6) + (2*x1*x5**2 + x1*(x5**2 - 3*x7**2) - 6*x3*x5*x7)*(2*x2*x6**2 + x2*(x6**2 - 3*x8**2) - 6*x4*x6*x8) + (6*x1*x5*x7 + 2*x3*x5**2 + x3*(x5**2 - 3*x7**2))*(6*x2*x6*x8 + 2*x4*x6**2 + x4*(x6**2 - 3*x8**2)))
    hess[4][6] =  2*(2*x1*(-sig_d + 2*x1*x5*x7 + 2*x2*x6*x8 + x3*(x5**2 - x7**2) + x4*(x6**2 - x8**2)) + 2*x3*(sig_c - x1*(x5**2 - x7**2) - x2*(x6**2 - x8**2) + 2*x3*x5*x7 + 2*x4*x6*x8) + 6*(x1*x5 - x3*x7)*(-sig_f + x1*x7*(3*x5**2 - x7**2) + x2*x8*(3*x6**2 - x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2)) + 6*(x1*x7 + x3*x5)*(sig_e - x1*x5*(x5**2 - 3*x7**2) - x2*x6*(x6**2 - 3*x8**2) + x3*x7*(3*x5**2 - x7**2) + x4*x8*(3*x6**2 - x8**2)) - (2*x1*x5**2 + x1*(x5**2 - 3*x7**2) - 6*x3*x5*x7)*(6*x1*x5*x7 - 2*x3*x7**2 + x3*(3*x5**2 - x7**2)) - (2*x1*x7**2 - x1*(3*x5**2 - x7**2) + 6*x3*x5*x7)*(6*x1*x5*x7 + 2*x3*x5**2 + x3*(x5**2 - 3*x7**2)))
    hess[4][7] =  2*(-x1*x4 + x2*x3 - 4*(x1*x5 - x3*x7)*(x2*x8 + x4*x6) + 4*(x1*x7 + x3*x5)*(x2*x6 - x4*x8) - (2*x1*x5**2 + x1*(x5**2 - 3*x7**2) - 6*x3*x5*x7)*(6*x2*x6*x8 - 2*x4*x8**2 + x4*(3*x6**2 - x8**2)) - (2*x2*x8**2 - x2*(3*x6**2 - x8**2) + 6*x4*x6*x8)*(6*x1*x5*x7 + 2*x3*x5**2 + x3*(x5**2 - 3*x7**2)))
    hess[5][0] =  2*(x2*x5 + x4*x7 + 4*x5*x7*(x2*x8 + x4*x6) + x5*(x5**2 - 3*x7**2)*(2*x2*x6**2 + x2*(x6**2 - 3*x8**2) - 6*x4*x6*x8) + x7*(3*x5**2 - x7**2)*(6*x2*x6*x8 + 2*x4*x6**2 + x4*(x6**2 - 3*x8**2)) + 2*(x5**2 - x7**2)*(x2*x6 - x4*x8))
    hess[5][1] =  2*(-sig_a + x1*x5 + 2*x2*x6 - x3*x7 + 4*x6*x8*(x2*x8 + x4*x6) + 6*x6*x8*(-sig_f + x1*x7*(3*x5**2 - x7**2) + x2*x8*(3*x6**2 - x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2)) + x6*(x6**2 - 3*x8**2)*(2*x2*x6**2 + x2*(x6**2 - 3*x8**2) - 6*x4*x6*x8) - 2*x6*(sig_c - x1*(x5**2 - x7**2) - x2*(x6**2 - x8**2) + 2*x3*x5*x7 + 2*x4*x6*x8) + x8*(3*x6**2 - x8**2)*(6*x2*x6*x8 + 2*x4*x6**2 + x4*(x6**2 - 3*x8**2)) + 2*x8*(-sig_d + 2*x1*x5*x7 + 2*x2*x6*x8 + x3*(x5**2 - x7**2) + x4*(x6**2 - x8**2)) + 2*(x6**2 - x8**2)*(x2*x6 - x4*x8) - 3*(x6**2 - x8**2)*(sig_e - x1*x5*(x5**2 - 3*x7**2) - x2*x6*(x6**2 - 3*x8**2) + x3*x7*(3*x5**2 - x7**2) + x4*x8*(3*x6**2 - x8**2)))
    hess[5][2] =  2*(-x2*x7 + x4*x5 - 4*x5*x7*(x2*x6 - x4*x8) + x5*(x5**2 - 3*x7**2)*(6*x2*x6*x8 + 2*x4*x6**2 + x4*(x6**2 - 3*x8**2)) - x7*(3*x5**2 - x7**2)*(2*x2*x6**2 + x2*(x6**2 - 3*x8**2) - 6*x4*x6*x8) + 2*(x5**2 - x7**2)*(x2*x8 + x4*x6))
    hess[5][3] =  2*(-sig_b + x1*x7 + x3*x5 + 2*x4*x6 - 4*x6*x8*(x2*x6 - x4*x8) + 6*x6*x8*(sig_e - x1*x5*(x5**2 - 3*x7**2) - x2*x6*(x6**2 - 3*x8**2) + x3*x7*(3*x5**2 - x7**2) + x4*x8*(3*x6**2 - x8**2)) + x6*(x6**2 - 3*x8**2)*(6*x2*x6*x8 + 2*x4*x6**2 + x4*(x6**2 - 3*x8**2)) + 2*x6*(-sig_d + 2*x1*x5*x7 + 2*x2*x6*x8 + x3*(x5**2 - x7**2) + x4*(x6**2 - x8**2)) - x8*(3*x6**2 - x8**2)*(2*x2*x6**2 + x2*(x6**2 - 3*x8**2) - 6*x4*x6*x8) + 2*x8*(sig_c - x1*(x5**2 - x7**2) - x2*(x6**2 - x8**2) + 2*x3*x5*x7 + 2*x4*x6*x8) + 2*(x6**2 - x8**2)*(x2*x8 + x4*x6) + 3*(x6**2 - x8**2)*(-sig_f + x1*x7*(3*x5**2 - x7**2) + x2*x8*(3*x6**2 - x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2)))
    hess[5][4] =  2*(x1*x2 + x3*x4 + 4*(x1*x5 - x3*x7)*(x2*x6 - x4*x8) + 4*(x1*x7 + x3*x5)*(x2*x8 + x4*x6) + (2*x1*x5**2 + x1*(x5**2 - 3*x7**2) - 6*x3*x5*x7)*(2*x2*x6**2 + x2*(x6**2 - 3*x8**2) - 6*x4*x6*x8) + (6*x1*x5*x7 + 2*x3*x5**2 + x3*(x5**2 - 3*x7**2))*(6*x2*x6*x8 + 2*x4*x6**2 + x4*(x6**2 - 3*x8**2)))
    hess[5][5] =  2*(x2**2 - 2*x2*(sig_c - x1*(x5**2 - x7**2) - x2*(x6**2 - x8**2) + 2*x3*x5*x7 + 2*x4*x6*x8) + x4**2 + 2*x4*(-sig_d + 2*x1*x5*x7 + 2*x2*x6*x8 + x3*(x5**2 - x7**2) + x4*(x6**2 - x8**2)) + 4*(x2*x6 - x4*x8)**2 - 6*(x2*x6 - x4*x8)*(sig_e - x1*x5*(x5**2 - 3*x7**2) - x2*x6*(x6**2 - 3*x8**2) + x3*x7*(3*x5**2 - x7**2) + x4*x8*(3*x6**2 - x8**2)) + 4*(x2*x8 + x4*x6)**2 + 6*(x2*x8 + x4*x6)*(-sig_f + x1*x7*(3*x5**2 - x7**2) + x2*x8*(3*x6**2 - x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2)) + (2*x2*x6**2 + x2*(x6**2 - 3*x8**2) - 6*x4*x6*x8)**2 + (6*x2*x6*x8 + 2*x4*x6**2 + x4*(x6**2 - 3*x8**2))**2)
    hess[5][6] =  2*(x1*x4 - x2*x3 + 4*(x1*x5 - x3*x7)*(x2*x8 + x4*x6) - 4*(x1*x7 + x3*x5)*(x2*x6 - x4*x8) - (2*x1*x7**2 - x1*(3*x5**2 - x7**2) + 6*x3*x5*x7)*(6*x2*x6*x8 + 2*x4*x6**2 + x4*(x6**2 - 3*x8**2)) - (2*x2*x6**2 + x2*(x6**2 - 3*x8**2) - 6*x4*x6*x8)*(6*x1*x5*x7 - 2*x3*x7**2 + x3*(3*x5**2 - x7**2)))
    hess[5][7] =  2*(2*x2*(-sig_d + 2*x1*x5*x7 + 2*x2*x6*x8 + x3*(x5**2 - x7**2) + x4*(x6**2 - x8**2)) + 2*x4*(sig_c - x1*(x5**2 - x7**2) - x2*(x6**2 - x8**2) + 2*x3*x5*x7 + 2*x4*x6*x8) + 6*(x2*x6 - x4*x8)*(-sig_f + x1*x7*(3*x5**2 - x7**2) + x2*x8*(3*x6**2 - x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2)) + 6*(x2*x8 + x4*x6)*(sig_e - x1*x5*(x5**2 - 3*x7**2) - x2*x6*(x6**2 - 3*x8**2) + x3*x7*(3*x5**2 - x7**2) + x4*x8*(3*x6**2 - x8**2)) - (2*x2*x6**2 + x2*(x6**2 - 3*x8**2) - 6*x4*x6*x8)*(6*x2*x6*x8 - 2*x4*x8**2 + x4*(3*x6**2 - x8**2)) - (2*x2*x8**2 - x2*(3*x6**2 - x8**2) + 6*x4*x6*x8)*(6*x2*x6*x8 + 2*x4*x6**2 + x4*(x6**2 - 3*x8**2)))
    hess[6][0] =  2*(-sig_b + 2*x1*x7 + x2*x8 + x4*x6 + 4*x5*x7*(x1*x5 - x3*x7) + 6*x5*x7*(sig_e - x1*x5*(x5**2 - 3*x7**2) - x2*x6*(x6**2 - 3*x8**2) + x3*x7*(3*x5**2 - x7**2) + x4*x8*(3*x6**2 - x8**2)) - x5*(x5**2 - 3*x7**2)*(6*x1*x5*x7 - 2*x3*x7**2 + x3*(3*x5**2 - x7**2)) + 2*x5*(-sig_d + 2*x1*x5*x7 + 2*x2*x6*x8 + x3*(x5**2 - x7**2) + x4*(x6**2 - x8**2)) - x7*(3*x5**2 - x7**2)*(2*x1*x7**2 - x1*(3*x5**2 - x7**2) + 6*x3*x5*x7) + 2*x7*(sig_c - x1*(x5**2 - x7**2) - x2*(x6**2 - x8**2) + 2*x3*x5*x7 + 2*x4*x6*x8) - 2*(x5**2 - x7**2)*(x1*x7 + x3*x5) + 3*(x5**2 - x7**2)*(-sig_f + x1*x7*(3*x5**2 - x7**2) + x2*x8*(3*x6**2 - x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2)))
    hess[6][1] =  2*(x1*x8 - x3*x6 + 4*x6*x8*(x1*x5 - x3*x7) - x6*(x6**2 - 3*x8**2)*(6*x1*x5*x7 - 2*x3*x7**2 + x3*(3*x5**2 - x7**2)) - x8*(3*x6**2 - x8**2)*(2*x1*x7**2 - x1*(3*x5**2 - x7**2) + 6*x3*x5*x7) - 2*(x6**2 - x8**2)*(x1*x7 + x3*x5))
    hess[6][2] =  2*(sig_a - x2*x6 + 2*x3*x7 + x4*x8 + 4*x5*x7*(x1*x7 + x3*x5) - 6*x5*x7*(-sig_f + x1*x7*(3*x5**2 - x7**2) + x2*x8*(3*x6**2 - x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2)) - x5*(x5**2 - 3*x7**2)*(2*x1*x7**2 - x1*(3*x5**2 - x7**2) + 6*x3*x5*x7) + 2*x5*(sig_c - x1*(x5**2 - x7**2) - x2*(x6**2 - x8**2) + 2*x3*x5*x7 + 2*x4*x6*x8) + x7*(3*x5**2 - x7**2)*(6*x1*x5*x7 - 2*x3*x7**2 + x3*(3*x5**2 - x7**2)) - 2*x7*(-sig_d + 2*x1*x5*x7 + 2*x2*x6*x8 + x3*(x5**2 - x7**2) + x4*(x6**2 - x8**2)) + 2*(x5**2 - x7**2)*(x1*x5 - x3*x7) + 3*(x5**2 - x7**2)*(sig_e - x1*x5*(x5**2 - 3*x7**2) - x2*x6*(x6**2 - 3*x8**2) + x3*x7*(3*x5**2 - x7**2) + x4*x8*(3*x6**2 - x8**2)))
    hess[6][3] =  2*(x1*x6 + x3*x8 + 4*x6*x8*(x1*x7 + x3*x5) - x6*(x6**2 - 3*x8**2)*(2*x1*x7**2 - x1*(3*x5**2 - x7**2) + 6*x3*x5*x7) + x8*(3*x6**2 - x8**2)*(6*x1*x5*x7 - 2*x3*x7**2 + x3*(3*x5**2 - x7**2)) + 2*(x6**2 - x8**2)*(x1*x5 - x3*x7))
    hess[6][4] =  2*(2*x1*(-sig_d + 2*x1*x5*x7 + 2*x2*x6*x8 + x3*(x5**2 - x7**2) + x4*(x6**2 - x8**2)) + 2*x3*(sig_c - x1*(x5**2 - x7**2) - x2*(x6**2 - x8**2) + 2*x3*x5*x7 + 2*x4*x6*x8) + 6*(x1*x5 - x3*x7)*(-sig_f + x1*x7*(3*x5**2 - x7**2) + x2*x8*(3*x6**2 - x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2)) + 6*(x1*x7 + x3*x5)*(sig_e - x1*x5*(x5**2 - 3*x7**2) - x2*x6*(x6**2 - 3*x8**2) + x3*x7*(3*x5**2 - x7**2) + x4*x8*(3*x6**2 - x8**2)) - (2*x1*x5**2 + x1*(x5**2 - 3*x7**2) - 6*x3*x5*x7)*(6*x1*x5*x7 - 2*x3*x7**2 + x3*(3*x5**2 - x7**2)) - (2*x1*x7**2 - x1*(3*x5**2 - x7**2) + 6*x3*x5*x7)*(6*x1*x5*x7 + 2*x3*x5**2 + x3*(x5**2 - 3*x7**2)))
    hess[6][5] =  2*(x1*x4 - x2*x3 + 4*(x1*x5 - x3*x7)*(x2*x8 + x4*x6) - 4*(x1*x7 + x3*x5)*(x2*x6 - x4*x8) - (2*x1*x7**2 - x1*(3*x5**2 - x7**2) + 6*x3*x5*x7)*(6*x2*x6*x8 + 2*x4*x6**2 + x4*(x6**2 - 3*x8**2)) - (2*x2*x6**2 + x2*(x6**2 - 3*x8**2) - 6*x4*x6*x8)*(6*x1*x5*x7 - 2*x3*x7**2 + x3*(3*x5**2 - x7**2)))
    hess[6][6] =  2*(x1**2 + 2*x1*(sig_c - x1*(x5**2 - x7**2) - x2*(x6**2 - x8**2) + 2*x3*x5*x7 + 2*x4*x6*x8) + x3**2 - 2*x3*(-sig_d + 2*x1*x5*x7 + 2*x2*x6*x8 + x3*(x5**2 - x7**2) + x4*(x6**2 - x8**2)) + 4*(x1*x5 - x3*x7)**2 + 6*(x1*x5 - x3*x7)*(sig_e - x1*x5*(x5**2 - 3*x7**2) - x2*x6*(x6**2 - 3*x8**2) + x3*x7*(3*x5**2 - x7**2) + x4*x8*(3*x6**2 - x8**2)) + 4*(x1*x7 + x3*x5)**2 - 6*(x1*x7 + x3*x5)*(-sig_f + x1*x7*(3*x5**2 - x7**2) + x2*x8*(3*x6**2 - x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2)) + (2*x1*x7**2 - x1*(3*x5**2 - x7**2) + 6*x3*x5*x7)**2 + (6*x1*x5*x7 - 2*x3*x7**2 + x3*(3*x5**2 - x7**2))**2)
    hess[6][7] =  2*(x1*x2 + x3*x4 + 4*(x1*x5 - x3*x7)*(x2*x6 - x4*x8) + 4*(x1*x7 + x3*x5)*(x2*x8 + x4*x6) + (2*x1*x7**2 - x1*(3*x5**2 - x7**2) + 6*x3*x5*x7)*(2*x2*x8**2 - x2*(3*x6**2 - x8**2) + 6*x4*x6*x8) + (6*x1*x5*x7 - 2*x3*x7**2 + x3*(3*x5**2 - x7**2))*(6*x2*x6*x8 - 2*x4*x8**2 + x4*(3*x6**2 - x8**2)))
    hess[7][0] =  2*(x2*x7 - x4*x5 + 4*x5*x7*(x2*x6 - x4*x8) - x5*(x5**2 - 3*x7**2)*(6*x2*x6*x8 - 2*x4*x8**2 + x4*(3*x6**2 - x8**2)) - x7*(3*x5**2 - x7**2)*(2*x2*x8**2 - x2*(3*x6**2 - x8**2) + 6*x4*x6*x8) - 2*(x5**2 - x7**2)*(x2*x8 + x4*x6))
    hess[7][1] =  2*(-sig_b + x1*x7 + 2*x2*x8 + x3*x5 + 4*x6*x8*(x2*x6 - x4*x8) + 6*x6*x8*(sig_e - x1*x5*(x5**2 - 3*x7**2) - x2*x6*(x6**2 - 3*x8**2) + x3*x7*(3*x5**2 - x7**2) + x4*x8*(3*x6**2 - x8**2)) - x6*(x6**2 - 3*x8**2)*(6*x2*x6*x8 - 2*x4*x8**2 + x4*(3*x6**2 - x8**2)) + 2*x6*(-sig_d + 2*x1*x5*x7 + 2*x2*x6*x8 + x3*(x5**2 - x7**2) + x4*(x6**2 - x8**2)) - x8*(3*x6**2 - x8**2)*(2*x2*x8**2 - x2*(3*x6**2 - x8**2) + 6*x4*x6*x8) + 2*x8*(sig_c - x1*(x5**2 - x7**2) - x2*(x6**2 - x8**2) + 2*x3*x5*x7 + 2*x4*x6*x8) - 2*(x6**2 - x8**2)*(x2*x8 + x4*x6) + 3*(x6**2 - x8**2)*(-sig_f + x1*x7*(3*x5**2 - x7**2) + x2*x8*(3*x6**2 - x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2)))
    hess[7][2] =  2*(x2*x5 + x4*x7 + 4*x5*x7*(x2*x8 + x4*x6) - x5*(x5**2 - 3*x7**2)*(2*x2*x8**2 - x2*(3*x6**2 - x8**2) + 6*x4*x6*x8) + x7*(3*x5**2 - x7**2)*(6*x2*x6*x8 - 2*x4*x8**2 + x4*(3*x6**2 - x8**2)) + 2*(x5**2 - x7**2)*(x2*x6 - x4*x8))
    hess[7][3] =  2*(sig_a - x1*x5 + x3*x7 + 2*x4*x8 + 4*x6*x8*(x2*x8 + x4*x6) - 6*x6*x8*(-sig_f + x1*x7*(3*x5**2 - x7**2) + x2*x8*(3*x6**2 - x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2)) - x6*(x6**2 - 3*x8**2)*(2*x2*x8**2 - x2*(3*x6**2 - x8**2) + 6*x4*x6*x8) + 2*x6*(sig_c - x1*(x5**2 - x7**2) - x2*(x6**2 - x8**2) + 2*x3*x5*x7 + 2*x4*x6*x8) + x8*(3*x6**2 - x8**2)*(6*x2*x6*x8 - 2*x4*x8**2 + x4*(3*x6**2 - x8**2)) - 2*x8*(-sig_d + 2*x1*x5*x7 + 2*x2*x6*x8 + x3*(x5**2 - x7**2) + x4*(x6**2 - x8**2)) + 2*(x6**2 - x8**2)*(x2*x6 - x4*x8) + 3*(x6**2 - x8**2)*(sig_e - x1*x5*(x5**2 - 3*x7**2) - x2*x6*(x6**2 - 3*x8**2) + x3*x7*(3*x5**2 - x7**2) + x4*x8*(3*x6**2 - x8**2)))
    hess[7][4] =  2*(-x1*x4 + x2*x3 - 4*(x1*x5 - x3*x7)*(x2*x8 + x4*x6) + 4*(x1*x7 + x3*x5)*(x2*x6 - x4*x8) - (2*x1*x5**2 + x1*(x5**2 - 3*x7**2) - 6*x3*x5*x7)*(6*x2*x6*x8 - 2*x4*x8**2 + x4*(3*x6**2 - x8**2)) - (2*x2*x8**2 - x2*(3*x6**2 - x8**2) + 6*x4*x6*x8)*(6*x1*x5*x7 + 2*x3*x5**2 + x3*(x5**2 - 3*x7**2)))
    hess[7][5] =  2*(2*x2*(-sig_d + 2*x1*x5*x7 + 2*x2*x6*x8 + x3*(x5**2 - x7**2) + x4*(x6**2 - x8**2)) + 2*x4*(sig_c - x1*(x5**2 - x7**2) - x2*(x6**2 - x8**2) + 2*x3*x5*x7 + 2*x4*x6*x8) + 6*(x2*x6 - x4*x8)*(-sig_f + x1*x7*(3*x5**2 - x7**2) + x2*x8*(3*x6**2 - x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2)) + 6*(x2*x8 + x4*x6)*(sig_e - x1*x5*(x5**2 - 3*x7**2) - x2*x6*(x6**2 - 3*x8**2) + x3*x7*(3*x5**2 - x7**2) + x4*x8*(3*x6**2 - x8**2)) - (2*x2*x6**2 + x2*(x6**2 - 3*x8**2) - 6*x4*x6*x8)*(6*x2*x6*x8 - 2*x4*x8**2 + x4*(3*x6**2 - x8**2)) - (2*x2*x8**2 - x2*(3*x6**2 - x8**2) + 6*x4*x6*x8)*(6*x2*x6*x8 + 2*x4*x6**2 + x4*(x6**2 - 3*x8**2)))
    hess[7][6] =  2*(x1*x2 + x3*x4 + 4*(x1*x5 - x3*x7)*(x2*x6 - x4*x8) + 4*(x1*x7 + x3*x5)*(x2*x8 + x4*x6) + (2*x1*x7**2 - x1*(3*x5**2 - x7**2) + 6*x3*x5*x7)*(2*x2*x8**2 - x2*(3*x6**2 - x8**2) + 6*x4*x6*x8) + (6*x1*x5*x7 - 2*x3*x7**2 + x3*(3*x5**2 - x7**2))*(6*x2*x6*x8 - 2*x4*x8**2 + x4*(3*x6**2 - x8**2)))
    hess[7][7] =  2*(x2**2 + 2*x2*(sig_c - x1*(x5**2 - x7**2) - x2*(x6**2 - x8**2) + 2*x3*x5*x7 + 2*x4*x6*x8) + x4**2 - 2*x4*(-sig_d + 2*x1*x5*x7 + 2*x2*x6*x8 + x3*(x5**2 - x7**2) + x4*(x6**2 - x8**2)) + 4*(x2*x6 - x4*x8)**2 + 6*(x2*x6 - x4*x8)*(sig_e - x1*x5*(x5**2 - 3*x7**2) - x2*x6*(x6**2 - 3*x8**2) + x3*x7*(3*x5**2 - x7**2) + x4*x8*(3*x6**2 - x8**2)) + 4*(x2*x8 + x4*x6)**2 - 6*(x2*x8 + x4*x6)*(-sig_f + x1*x7*(3*x5**2 - x7**2) + x2*x8*(3*x6**2 - x8**2) + x3*x5*(x5**2 - 3*x7**2) + x4*x6*(x6**2 - 3*x8**2)) + (2*x2*x8**2 - x2*(3*x6**2 - x8**2) + 6*x4*x6*x8)**2 + (6*x2*x6*x8 - 2*x4*x8**2 + x4*(3*x6**2 - x8**2))**2)
    return hess
