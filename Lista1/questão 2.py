
import matplotlib.pyplot as plt
import numpy as np
from math import *
from cmath import *


#-du²/d²x = 1 - u / - ε
#u(0) = 0
#u(1) = 0


#axes definitions
ax = plt.subplot(1,1,1)

ax.grid()

#Coeficientes

c1 
c2 
E  #Constante e deve ser diferente de zero

c2 = (exp((-1)/sqtr(E)) - 1) / (exp(1/sqtr(E)) - exp((-1)/sqtr(E))) 

c1 = - 1 - c2

#Definição da função
def F(x): 
    return c1 * exp(- x / sqtr (E)) +  c2 * exp(x / sqtr (E)) + 1


#Solicitar o valor de E
#testes com E = 0.1, 0.01, 0.001, 0.0001
E = 0.1


#analitica
#Expansão com formula de Taylor até a 2 ordem

xx = np.linspace (0, 1)
yy = np.zeros_like(xx) #retorna uma matriz de zeros
for i,x in enumerate(xx):
    yy[i] = F(x)
    

#Definição numerica
def un(N):
#malha
   a = 0
   b = 1
   h = (b-a)/(N-1)
   xp = np.linspace(a,b,N) #começa em 0 e termina em 1

   A = np.zeros((N,N)) #retorna uma matriz de zeros
   b = np.zeros(N) #retorna uma matriz de zeros

   A[0,0] = 1
   b[0] = 0
   for i in np.arange(1,N-1):
       A[i,i-1] = 1
       A[i,i] = -2
       A[i,i+1] = 1
       b[i] = #Formula com o polinomio de Taylor
   A[N-1,N-1] = 1
   b[N-1] = 0

   u = np.linalg.solve(A,b)
   return (xp,u)

xx, F = un(11)
ax.plot(xx, F, 'ro', markersize=4)
