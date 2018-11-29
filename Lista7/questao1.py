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
figure(num=None, figsize=(8, 6), dpi=72, facecolor='w', edgecolor='k')


###### Testes #######

#funcao do slide

ElementosFinitos().funcaoBase()