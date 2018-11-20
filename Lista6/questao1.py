# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 16:39:54 2018

@author: thiagoalmeida
"""

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.pyplot import figure
import numpy as np

from MetodosIntegracao import NewtonCotes
from MetodosIntegracao import Repetidos


#define o tamanho dos graficos
figure(num=None, figsize=(8, 6), dpi=72, facecolor='w', edgecolor='k')


###### Testes #######

#funcao do slide
def f(x):
    return np.exp(x) * np.cos(x)

#intervalo
a = 0.0
b = 1.2
num_divisoes = 120
num_divisoes2 = 2*60
num_divisoes3 = 3*40

print("Valor exato", 1.648774427347253)
print("NewtonCotes().Retangulo", NewtonCotes().Retangulo(a, b, f))
print("NewtonCotes().PontoMedio", NewtonCotes().PontoMedio(a, b, f))
print("NewtonCotes().Trapezio", NewtonCotes().Trapezio(a, b, f))
print("NewtonCotes().Simpsom13", NewtonCotes().Simpsom13(a, b, f))
print("NewtonCotes().Simpsom38", NewtonCotes().Simpsom38(a, b, f))
print("--- repetidos ---")
print("Repetidos().Retangulo", Repetidos().Retangulo(a, b, num_divisoes, f))
print("Repetidos().PontoMedio", Repetidos().PontoMedio(a, b, num_divisoes, f))
print("Repetidos().Trapezio", Repetidos().Trapezio(a, b, num_divisoes, f))
print("Repetidos().Simpsom13", Repetidos().Simpsom13(a, b, num_divisoes2, f))
print("Repetidos().Simpsom38", Repetidos().Simpsom38(a, b, num_divisoes3, f))

###### Funcoes da lista #######


def u(c1, c2, E, x):
    t = x / np.sqrt(E)
    return c1 * np.exp(-t) + c2 * np.exp(t) + 1.0
    
    
    
###### Definicoes #######
    
# Funcao u(x)    
E = 0.001
t = 1.0 / np.sqrt(E)
c1 = -1.0
c2 = (np.exp(-t) - 1.0) / (np.exp(t) - np.exp(-t))

#intervalo de integracao
a = 0.0
b = 1.0