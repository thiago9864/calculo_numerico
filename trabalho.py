# -*- coding: utf-8 -*-

from math import *
from cmath import *
import matplotlib.pyplot as plt

#import pandas as pd
#import queue
#Definir as variaveis

global T0 #Tempo inicial
global Tm #Tempo final

global K #constante da formula
K = 0.035871952 

global C #constante para salvar parte da formula

T0 = 50
Tm = 0.01

C = (T0-Tm) 

global t
global ti
global n

t = 5
ti = 0.01

#Listas
euler_explicito = [] #lista para salvar os dados da primeira função
euler_implicito = [] #lista para salvar os dados da segunda função
intervalo = []

#Definição das funções

def solucaoInicial(): #Solução inicial do problema
    return (C*exp((-1)*K*t)+Tm)

def funcao_original(X,t): #Função Original do Problema
    return ((-1)*K*C)

def _h(): #variação de um intervalo, no caso o tempo
    return (t-ti)/n

    #Definição do Metodo de Euler Explícito
def Euler_Explicito():
    y = solucaoInicial()
    x = ti
    h = _h()
    
    i = 0
    while i < n:
        y = y + h * funcao_original(x,y)
        x = x + h
        i = i + 1
        euler_explicito.append(y)
        intervalo.append(x)
    
    return y

n = 64#int(input("O número de subdivisões: ")) #solicitar o numero de divisões
print("teste")
Euler_Explicito()
plt.plot(intervalo, euler_explicito)
plt.ylabel('some numbers')
plt.show()

#Definição do Metodo de Euler Implícito
def euler_implicito():
    y_imp = solucaoInicial()
    x_imp = ti
    h_imp = h()
    
    y_p,x_p = 0,0
    
    i = 0
    while i < n:
        y_p = y_imp + h_imp*f(x_imp, y_imp)
        X_p = x_imp + h_imp
        y_imp = y_imp + (h_imp/2)*(f(x_imp, y_imp) + f(x_p,y_p))
        x_imp = x_imp + h_imp
        i = i + 1
        
    euler_implicito.append(y_imp)
    
    return y_imp