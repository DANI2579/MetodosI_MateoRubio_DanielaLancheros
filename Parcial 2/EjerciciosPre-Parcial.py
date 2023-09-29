################################################################################################################################################################
########################################################################Punto 17################################################################################
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import sympy as sym

e = np.e

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
        np.append(Roots,0.)
        
    Roots.sort()
    
    return Roots
    
def GetAllRootsGLag(n):
    x = sym.Symbol('x',real = True)
    X = np.linspace(0,n+((n-1)*np.sqrt(n)),100)
    pol = GetLaguerreRecursive(n,x)
    Poly = sym.lambdify([x],pol,'numpy')
    Dpoly = sym.lambdify([x],sym.diff(pol,x,1),'numpy')
    return GetRoots(Poly,Dpoly,X)

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

def GetPunto17():
    f = lambda x: (x**3)*(e**x)/(((e**x)-1))
    N = np.zeros(9)
    v = np.zeros(9)
    for i in range(2,11):
        N[i-2] = i
        v[i-2] = IntegrationLaguerre(i,f)/(np.pi**4/15)
    plt.scatter(N,v)
    plt.xlabel("n")
    plt.ylabel("Error Relativo")
    plt.title("Grafica de Error relativo para el metodo de Lagerre")
    plt.axhline(y = 1,color='r')
    plt.show()  
    return IntegrationLaguerre(3,f)

#print(GetPunto17())

################################################################################################################################################################
########################################################################Punto 18################################################################################

def GetHermite(n,x):
    return ((-1)**n)*(e**(x**2))*(sym.diff((e**(-(x**(2)))),x,n))


def GetAllRootsGHer(n):
    x = sym.Symbol('x',real = True)
    X = np.linspace(-np.sqrt((4*(n))+1),np.sqrt((4*(n))+1),10000)
    pol = GetHermite(n,x)
    Poly = sym.lambdify([x],pol,'numpy')
    Dpoly = sym.lambdify([x],sym.diff(pol,x,1),'numpy')
    return GetRoots(Poly,Dpoly,X)

def GetWeightsGHer(n):
    roots = GetAllRootsGHer(n)
    a = (2**(n-1))*np.math.factorial(n)*np.sqrt(np.pi)
    return a/((n**2)*((roots*GetHermite(n-1,roots))**2))


#n = 1
#print(GetAllRootsGHer(n))
#print(np.polynomial.hermite.hermgauss(n)[0])