from sympy import diff, Symbol, symbols, exp
x1,x2,x3,x4,x5 = symbols('x1 x2 x3 x4 x5')
yi, ti = symbols('yi ti')
ri = yi - (x1 + x2 * exp(-ti*x4) + x3 * exp(-ti*x5))

#result = open('res_dipole','w')
x_array = [x1,x2,x3,x4,x5]
#for i,x in enumerate(x_array):
    #print 'der[' + str(i) + '] = ',
    #print diff(f,x)
for i,x in enumerate(x_array):
    for j,y in enumerate(x_array):
        print 'hess[' + str(i) + '][' + str(j) + '] = ',
        print diff(ri,x,y)

