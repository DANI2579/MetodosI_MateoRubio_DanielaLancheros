import numpy as np
import matplotlib.pyplot as plt
#######Inciso a)

def diferencia_dividida(X,Y):
    matrix = np.zeros((X.shape[0],Y.shape[0]))
    matrix[:,0] = Y
    for j in range(1,matrix.shape[0]):
        for i in range(matrix.shape[1]-j):
            matrix[i][j] = (matrix[i+1][j-1] - matrix[i][j-1])/(X[i+j]-X[i])
    return matrix

def f(x): 
    return np.power(np.e,-x) - x

X = np.linspace(0,2,3)
Y = f(X)
n3 = diferencia_dividida(X,Y)[0][2]
n2 = diferencia_dividida(X,Y)[0][1]

a = n3
b = n2 - ((X[0] + X[1])*n3)
c = Y[0] - (X[0]*n2) + (X[0]*X[1]*n3)

result1 = (-2*c)/(b + np.sqrt((b**2) - 4*a*c))
result2 = (-2*c)/(b - np.sqrt((b**2) - 4*a*c))
print("Resltado inciso a) " +  str(result2))

#######################################################
##Cuando f es evaluada en -1 el valor de f es e + 1 que es positivo
# Cuando f es evaluada en 1 el valor de f es e^-1 -1 que es negativo
# Por tanto se concluye que entre -1 y 1 debe existir una raiz para la funcion dada
x0 = 1
x1 = -1
print("El resultado del inciso b) son x0 = {}, x1 = {}".format(x0,x1))
#######################################################
x2 = (x0+ x1)/2
print("x2 = " + str(x2)) 
#######################################################
##Dado el algoritmo el valor esta en x3
if f(x2)*f(x1) > 0:
    x3 = x2
else: 
    x3 = x2
#######################################################
def biseccion(f,a,b,itmax =100, preci = 1e-10):
    c = (a+b)/2          
    i = 0 
    while np.abs(f(c)) > preci and i<itmax:
        c = (a+b)/2
        if f(c)*f(a) > 0:
            a = c
        else: 
            b = c
    return c

print("La raiz de la funcion es: " + str(biseccion(f,x1,x0)))












