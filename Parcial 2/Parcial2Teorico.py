import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import sympy as sym

#####Punto A ######################################################################

def GetLaguerreRecursive(n,x):
    '''
    Devuelve los polinomios de laguerre dado un x
    '''
    if n == 0:
        return 1
    elif n == 1: 
        return 1-x 
    else: 
        return ((((2*(n-1))+(1-x))*GetLaguerreRecursive(n-1,x))-((n-1)*GetLaguerreRecursive(n-2,x)))/(n)
    
x = sym.Symbol("x", Real = True)
print(sym.simplify(GetLaguerreRecursive(2,x)))
#####Punto B ############################################################
def GetNewton(f,df,xn,itmax=100,precision=1e-14):
    
    error = 1.
    it = 0
    
    while error >= precision and it < itmax:
        
        try:
            
            xn1 = xn - f(xn)/df(xn)
            
            error = np.abs(f(xn)/df(xn))
            
        except ZeroDivisionError:
            print('Zero Division')
            
        xn = xn1
        it += 1
        
    if it == itmax:
        return False
    else:
        return xn
    
def GetRoots(f,df,x,tolerancia = 13):
    
    Roots = np.array([])
    
    for i in x:
        
        root = GetNewton(f,df,i)
        
        if root != False:
            
            croot = np.round( root, tolerancia )
            
            if croot not in Roots:
                Roots = np.append(Roots, croot)
    if np.abs(f(0)) == 0. and not 0. in Roots:
        Roots = np.append(Roots,0.)
        
    Roots.sort()
    
    return Roots
    
def GetAllRootsGLag(n):
    x = sym.Symbol('x',real = True)
    X = np.linspace(0,n+((n-1)*np.sqrt(n)),100)
    pol = GetLaguerreRecursive(n,x)
    Poly = sym.lambdify([x],pol,'numpy')
    Dpoly = sym.lambdify([x],sym.diff(pol,x,1),'numpy')
    return GetRoots(Poly,Dpoly,X)



##### Punto B ######################################################################
Roots = GetAllRootsGLag(2)
r1 = Roots[0]
r2 = Roots[1]
print("El resultado del punto B) es " + str(Roots))

###### Punto C######################################################################
def GetWeightsGLag(n): 
    roots = GetAllRootsGLag(n)
    return (roots)/(((n+1)**2)*((GetLaguerreRecursive(n+1,roots))**2))

def IntegrationLaguerre(n,f):
    I = 0
    pesos = GetWeightsGLag(n)
    raices = GetAllRootsGLag(n)
    for i in range(n):
        I += pesos[i]*f(raices[i])
    return I

f = lambda x:((x-Roots[1])/(Roots[0]-Roots[1]))
f1 = lambda x:((x-Roots[0])/(Roots[1]-Roots[0]))
w1 = IntegrationLaguerre(2,f)
w2 = IntegrationLaguerre(2,f1)
print( "El resultado del punto C)) es"+ "Peso 1: {}, Peso 2: {}".format(round(w1,8),round(w2,8)))

###################### Punto D #############################################
f = lambda x: x**3
Resultado = w1*f(r1)
Resultado += w2*f(r2)
print("Dado que la regla de Gamma es (x-1)! y reemplazando x= 4, entonces el resultado de la integral es 6")
print("La respuesta del punto D) es " + str(Resultado))
    
