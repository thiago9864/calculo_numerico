
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



