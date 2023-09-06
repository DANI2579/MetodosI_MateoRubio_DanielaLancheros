import numpy as np
import matplotlib.pyplot as plt



def derivada_polinomio(f,x,h=0.01):
    return (1/(2*h))*(f(x)*(-3) + 4*f(x+h) - f(x+(2*h)))

def derivada_centrada(f,x,h = 0.01):
    if h != 0:
        return(f(x+h)-f(x-h))/(2*h)
    
def derivada_exacta(x):
    return (1/2)*(1/np.sqrt(np.tan(x)))*np.square(1/np.cos(x))

def root_tan(x): 
    return(np.sqrt(np.tan(x)))

x = np.arange(0.1,1.1,0.01)
df_dxE = derivada_exacta(x)
df_dxP = derivada_polinomio(root_tan,x)
df_dxC = derivada_centrada(root_tan,x)
plt.suptitle("Graficacion de distintos metodos de derivacion")
plt.plot(x,df_dxE, label = 'Derivada Exacta')
plt.scatter(x,df_dxP,color ='green', s = 10 ,label = 'Derivada Polinomio')
plt.scatter(x,df_dxC, color ='orange', s = 10,label = 'Derivada Centrada')
plt.legend()

plt.show()
plt.close()
plt.suptitle("Graficacion del error nodal de distintos metodos de derivacion")
plt.plot(x,np.abs(df_dxE- df_dxE),label = 'Derivada Exacta')
plt.plot(x,np.abs(df_dxE- df_dxP),label = 'Derivada Polinomica')
plt.plot(x,np.abs(df_dxE- df_dxC),label = 'Derivada Centrada')
plt.legend()
plt.show()

##Rta: Si, efectivamente ambas aproximaciones tienen un orden de precision dado por
# el mismo exponente 
