import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import sympy as sym

e = np.e
h = 6.626*(10**-34)
k = 1.3806*(10**-23)
c = 3*(10**8)
T = 5772
l0 = 100*(10**-9)
l1 = 400*(10**-9)


roots, weights = np.polynomial.legendre.leggauss(20)
roots1,weights1 = np.polynomial.laguerre.laggauss(20)
v1 = c/l1
v0 = c/l0
b = (h*v0)/(k*T)
a = (h*v1)/(k*T)
f = lambda x: (x**3)/((np.e**x)-1) 
f1 = lambda x: (x**3)*(e**x)/(((e**x)-1))

def integralab(f,a,b,n):
    roots, weights = np.polynomial.legendre.leggauss(20)
    integral=(b-a)/2
    I=0
    for i in range(n):
        I+=weights[i]*f(((roots[i]*(b-a))/2)+((b+a)/2))
    integral*=I
    return integral

def IntegrateDenom(f):
    I = 0
    for i in range(20):
        I += weights1[i]*f(roots1[i])
    return I

resultado = integralab(f,a,b,20)/IntegrateDenom(f1) 
print(resultado)

#####D) Esto se debe a que gracias a la capa de ozono y la atmosfera,
###### los rayos ultravioletas se dispersan haciendo que le porcentaje
###### que llegue a la tierra sea menor. 
