# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 11:40:34 2018

@author: Thiago, Renan
"""
import numpy as np
import math as m

num_passos_bis = 0
num_passos_fp = 0

### Funcao pra achar a raiz ###

def F(x):
    return x*3 + 4*x**2 - 10
    
### Funcoes de ponto fixo dadas
    
def f1(x):
    return x - x**3 + 4*x**2 + 10
    
def f2(x):
    return (1/2) * m.sqrt(10 - x**3)
    
def f3(x):
    return x.sqrt(10 / (4+x))
    
def f4(x):
    return x - ((x**3 + 4*x**2 - 10) / (3*x**2 + 8*x))