import numpy as np
import matplotlib.pyplot as plt

def newton_raphson(f,df_dx,xi,itmax = 100,preci = 1e-8):
    if df_dx(f,xi,1e-6) != 0:
        while itmax-1 >0 and np.abs(f(xi)/df_dx(f,xi,1e-6))>preci:
            xi = xi-(f(xi)/df_dx(f,xi,1e-6))
            itmax-=1
    else:
        print("Division por cero")
        return False
    return xi

def Todas_las_raices(f,x,df_dx, aprox=5):
    resultado = np.array([])
    for i in x:
        raiz = newton_raphson(f,df_dx,i)
        if raiz != False:
            raiz = np.round(raiz,aprox)
            if raiz not in resultado:
                resultado = np.append(resultado,raiz)
    return resultado 

def df_dx_centrada(f,x,h):
    if h != 0:
        return(f(x+h)-f(x-h))/(2*h)

def f(x): 
    return 3*(x**5) + 5*(x**4) - (x**3)

x = np.linspace(-3,1,1000)
Todas_las_raices(f,x,df_dx_centrada) 

