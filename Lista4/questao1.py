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

def gerarPontos1(a, b, n):
    xk = np.zeros((n,), dtype=np.float128)
    
    for k in range (0,n):
        xk[k] = a + ((b-a)*k)/n        
    
    return xk
    
def gerarPontos2(a, b, n):
    xk = np.zeros((n,), dtype=np.float128)
    
    for k in range (0,n):
        xk[k] = ((a+b) / 2.0) - ((b-a) / 2.0) * m.cos(float(k)/n * m.pi)
    
    return xk
        
##### Gerador de grafico #####

def gerarGrafico(tempo, solucao_aproximada, solucao_exata, metodo):
    
    #print("len tempo: ", len(tempo))
    #print("len solucao_aproximada: ", len(solucao_aproximada))
    #print("len solucao_exata: ", len(solucao_exata))
    
    
     #plota grafico da função
    plt.plot(
        tempo, solucao_exata, 'b--',
        tempo, solucao_aproximada, 'r-'    
        )
    plt.ylabel(u"Valor de u(h)") #esse 'u' antes da string é pra converter o texto pra unicode
    plt.xlabel(u"Valor de h, " + str(len(tempo)) + u" partições")
    
    #legendas do grafico
    se_line = mlines.Line2D([], [], color='blue', marker='', markersize=0, label=u'Sol Exata')
    ac_line = mlines.Line2D([], [], color='red', marker='', markersize=0, label=u'Met. ' + metodo)
    
    plt.legend(handles=[se_line, ac_line])
    
    plt.title(u"Metodo "+metodo+u" x Solução exata", )
    
    #plt.axis([0, 50, 0, 100])
    plt.show()   

#configura grafico
n = 100
a = -1.0
b = 1.0

pontos = gerarPontos2(a, b, n)
sol_exata = np.zeros((n,), dtype=np.float128)
sol_aprox = np.zeros((n,), dtype=np.float128)
    
### teste ###
    
def f(x):
    return m.cos(x)

x = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7], np.float64)
#x = np.array([0.2, 0.3, 0.4], np.float64)

#############

print("")
print("------")
print("")

#print("Metodo de Gauss (direto)")
inicio = Utils().getTime()
Newton().interpolacaoIterativa(x, f, 3);
#print("Res Newton: " + repr(res));
fim  = Utils().getTime()




for k in range (0,n):
    sol_exata[k] = F(pontos[k])

gerarGrafico(pontos, sol_aprox, sol_exata, "Teste")