import matplotlib.pyplot as plt
import numpy as np
from math import *
from cmath import *

#Array
arr_problema_valor_contorno = []

#Coeficientes
c1 
c2 
E  #Constante e deve ser diferente de zero
Nel #numeros de elementos da malha

c1 = - 1 - c2

c2 = (exp((-1)/sqtr(E)) - 1) / (exp(1/sqtr(E)) - exp((-1)/sqtr(E))) 

valoresIniciais = 0 # u(0) = u(1) = 0

#Definindo as funções
def solucao_exata(): #solução exata dada
    return c1 * exp( ((-1) * x) / sqtr(E)) +   c2 * exp( x / sqtr(E)) + 1

def solucao_original(u): #derivada segunda original
    return (1 - u) / ((-1) * E)

def h(a,b): #Partição do Domínio em partes iguais
    h = (b - a)/ Nel

#Solicitar o valor de E
#testes com E = 0.1, 0.01, 0.001, 0.0001
E = 0.1
