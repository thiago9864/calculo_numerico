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

print("Valor exato", 1.648774427)
print("NewtonCotes().Retangulo", NewtonCotes().Retangulo(a, b, f))
print("NewtonCotes().PontoMedio", NewtonCotes().PontoMedio(a, b, f))
print("NewtonCotes().Trapezio", NewtonCotes().Trapezio(a, b, f))
print("NewtonCotes().Simpsom13", NewtonCotes().Simpsom13(a, b, f))
print("NewtonCotes().Simpsom38", NewtonCotes().Simpsom38(a, b, f))