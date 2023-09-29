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
    x = sym.Symbol('x',real = True)
    roots = GetAllRootsGHer(n)
    pol = GetHermite(n-1,x)
    Poly = sym.lambdify([x],pol,'numpy')
    a = (2**(n-1))*np.math.factorial(n)*np.sqrt(np.pi)
    return a/((n**2)*((Poly(roots))**2))

def IntegrateHermite(f,n):
    I=0
    pesos = GetWeightsGHer(n)
    functions = f(GetAllRootsGHer(n))
    for i in range(n):
        I+= pesos[i]*functions[i]
    return I
def GetPunto18():
    print(GetWeightsGHer(20))
    print(GetAllRootsGHer(20))
    f = lambda x: 4*(x**2)*(x**2)
    return  ((1/(np.sqrt(2)*(pow(np.pi,1/4))))**2)*IntegrateHermite(f,5) 

#print(GetPunto18())

################################################################################################################################################################
########################################################################Punto 19################################################################################

def GetLegendre(n,x):
    return  sym.diff( (x**2 - 1)**n,x,n )/(2**n*np.math.factorial(n))

def GetAllRootsGLeg(n):
    x = sym.Symbol('x',real = True)
    X = np.linspace(-1,1,100)
    pol = GetLegendre(n,x)
    Poly = sym.lambdify([x],pol,'numpy')
    Dpoly = sym.lambdify([x],sym.diff(pol,x,1),'numpy')
    return GetRoots(Poly,Dpoly,X)

def GetWeightsGLeg(n):
    x = sym.Symbol('x',real = True)
    Roots = GetAllRootsGLeg(n)
    pol = GetLegendre(n,x)
    Dpoly = sym.lambdify([x],sym.diff(pol,x,1),'numpy')
    return 2/((1-(np.square(Roots))) * (np.square(Dpoly(Roots))))

def IntegrateLegendre(f,n):
    I = 0
    pesos = GetWeightsGLeg(n)
    raices = GetAllRootsGLeg(n)
    for i in range(n): 
        I += pesos[i] * f(raices[i])
    return I

def GetPunto19():
    dT = 10**-4  
    i =0
    centinela = True 
    f = lambda x,T,delta : np.tanh(np.sqrt(np.square(x) + np.square(delta))*(300/(2*T)))/(np.sqrt(np.square(x) + np.square(delta)))
    while i <= 1 and centinela:
        fx = lambda x: f(x,12.1331,0)
        I = IntegrateLegendre(fx,34)
        if(np.abs(I-(1/(0.3))) <dT):
            centinela = False
        else:
            dT*=dT
        i+=1
    return I
#f = lambda x: x
#print(IntegrateLegendre(f,50))
print(GetPunto19())

#print(np.size(np.polynomial.legendre.leggauss(50)[0]))
#print(np.size(GetAllRootsGLeg(50)))

