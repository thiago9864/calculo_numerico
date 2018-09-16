# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import math as m
from array import array

### Funcao pra achar a raiz ###

def F(x):
    #return (x/2.0)**2 - m.sin(x)
    return x**2 - 3.0 * x + 2.0
    #return m.cos(x)

### Constantes ###

a = -5
b = 5
E = 0.001

## checa intervalos e procura menor com uma raiz ##
## se encontrar, coloca no array e parte pro proximo ##

intervalo = b - a

multiplicador = 2.0
np = intervalo * multiplicador
dx = intervalo / np

arr_intervalos = []

for i in range(1, int(np)):
    pa = ((i - 1) * dx) + a
    p = (i * dx) + a

    ya = F(pa)
    y = F(p)

    #print ("fx: ", ya, y)

    if((ya >= 0 and y < 0) or (ya < 0 and y >= 0)):
        arr_intervalos.append([pa, p])

def bissecao(_a, _b):

    e = _b - _a
    mx = e / 2.0
    _x = mx + _a
    iteracoes = 0

    #checa se a raiz esta nos intervalos
    if(F(_a) == 0.0):
        return _a
    if(F(_b) == 0.0):
        return _b

    while(e > E and iteracoes < 200):
        iteracoes = iteracoes + 1

        #acha meio
        e = _b - _a
        mx = e / 2.0
        _x = mx + _a

        #determina valores
        fa = F(_a)
        fb = F(_b)
        fx = F(_x)
        
        #se a raiz estiver entre um dos intervalos formados, ajusta intervalo (_a, _b)
        if((fa >= 0 and fx < 0) or (fa < 0 and fx >= 0)):
            #intervalo a-x
            _b = _x
        else:
            #intervalo x-b
            _a = _x

    #atingiu condicao de parada
    return _x

def FalsaPosicao(_a,_b):
    
    e = E #epsilon - precisão fornecida
    
    #Checa a condição de existencia
    if (F(_a) * F(_b) < 0):
        
        while (abs(F(_x)) > e and abs(_b - _a) > e):
            
            if(F(_a) * F(_b) < 0):
                _a = _x
            else:
                _b = _x
                
            _x = (_a * F(_b) + _b * F(_a)) / (F(_b) - F(_a))
        
        return _x

#print(arr_intervalos)

if(len(arr_intervalos) == 0):
    print("Nenhuma raiz encontrada")
else:
    for i in range(len(arr_intervalos)):
        intervalo = arr_intervalos[i]
        print("raiz #" + str(i))
        print ("pesquisar no intervalo (" + str(intervalo[0]) + ", " + str(intervalo[1]) + ")")
        raiz = bissecao(intervalo[0], intervalo[1])
        print("raiz encontrada: " + str(raiz) + "\n")

