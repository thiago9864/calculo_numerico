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

a = 0
b = 1
N = 4
ordem_matriz = 2

x = np.zeros((N,), dtype=np.float128)

h = (b-a)/N
d = 0.0

for i in range(N):
    d += h
    x[i] = d

print("h:", h)
print("x:", x)

#faz a matriz pequena
K = np.zeros((ordem_matriz,ordem_matriz), dtype=np.float128)
for i in range(ordem_matriz):
    for j in range(ordem_matriz):
        K[i][j] = ElementosFinitos().K(0, 0, h, N, x) * (1.0 / h)

print("K:", K)