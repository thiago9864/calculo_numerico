# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 11:12:42 2018

@author: thiagoalmeida
"""


import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.pyplot import figure
import numpy as np

from MetodosIntegracao import QuadraturaGauss
from Erros import Erro
from ElementosFinitos import ElementosFinitos


#define o tamanho dos graficos
#figure(num=None, figsize=(8, 6), dpi=72, facecolor='w', edgecolor='k')


###### Testes #######

#funcao do slide

a = 0.0
b = 1.0
elementos = 4
grau_polinomio = 1


pontos_elementos = elementos + 1
pontos_polinomio = grau_polinomio + 1

x = np.zeros((pontos_elementos,), dtype=np.float128)

h = (b-a)/elementos
d = a

for i in range(pontos_elementos):
    x[i] = d
    d += h

print("h:", h)
print("x:", x)

#faz a matriz pequena
K = np.zeros((pontos_polinomio,pontos_polinomio), dtype=np.float128)
for i in range(pontos_polinomio):
    for j in range(pontos_polinomio):
        K[i][j] = ElementosFinitos().K(i, j, h, pontos_polinomio)
        
#faz o termo fonte pequeno
F = np.zeros((pontos_polinomio,), dtype=np.float128)
for i in range(pontos_polinomio):
        F[i] = ElementosFinitos().F(i, h, pontos_polinomio)

print("K:", K)
print("F:", F)