#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 21:10:01 2018

@author: Thiago Almeida
"""

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import math as m

#imports locais
from Utils import Utils
#from interpoladores.Lagrange import Lagrange
from interpoladores.Newton import Newton


#funcao exata
def F(x):
    return 1.0 / (1.0 + 25.0 * x**2)


##### Gerador de grafico #####

def f(x):
    return m.cos(x);

x = np.array([0.2, 0.3, 0.4], np.float64)

print("")
print("------")
print("")

#print("Metodo de Gauss (direto)")
inicio = Utils().getTime()
res = Newton().executar(x, f, 3);
print("Res Newton: " + repr(res));
fim  = Utils().getTime()

