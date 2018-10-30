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
from interpoladores.Lagrange import Lagrange
from interpoladores.Linear import Linear
from interpoladores.Newton import Newton


#funcao exata
def F(x):
    #return m.cos(x)
    return 1.0 / (1.0 + 25.0 * x**2)

#gera pontos para letra A
def gerarPontos_a(a, b, n):
    xk = np.zeros((n+1,), dtype=np.float128)
    
    for k in range (0,n+1):
        xk[k] = a + ((b-a)*k)/n        
    
    return xk
    
#gera pontos para letra B
def gerarPontos_b(a, b, n):
    xk = np.zeros((n+1,), dtype=np.float128)
    
    for k in range (0,n+1):
        xk[k] = ((a+b) / 2.0) - ((b-a) / 2.0) * m.cos(float(k)/n * m.pi)
    
    return xk

#gera pontos para letra C
def gerarPontos_c(a, b, n):
    xk = np.zeros((n+1,), dtype=np.float128)
    
    for k in range (0,n+1):
        xk[k] = -1.0 + (2.0 * k) / float(n)
    
    return xk
        
##### Geradores de graficos #####

def gerarGrafico(tempo, solucao_aproximada, solucao_exata, metodo, ordem):
    
    #print("len tempo: ", len(tempo))
    #print("len solucao_aproximada: ", len(solucao_aproximada))
    #print("len solucao_exata: ", len(solucao_exata))
    
    
     #plota grafico da função
    plt.plot(
        tempo, solucao_exata, 'b--',
        tempo, solucao_aproximada, 'r-'    
        )
    #plt.ylabel(u"Valor de f(x), P(x)") #esse 'u' antes da string é pra converter o texto pra unicode
    #plt.xlabel(u"Valor de x, " + str(len(tempo)) + u" partições")
    plt.xlabel(str(len(tempo)) + u" partições, N = " + str(ordem))
    
    #legendas do grafico
    se_line = mlines.Line2D([], [], color='blue', marker='', markersize=0, label=u'Sol Exata')
    ac_line = mlines.Line2D([], [], color='red', marker='', markersize=0, label=u'Met. ' + metodo)
    
    plt.legend(handles=[se_line, ac_line])
    
    plt.title(u"Metodo "+metodo+u" x Solução exata", )
    
    #plt.axis([0, 50, 0, 100])
    plt.show()  
    
def gerarGraficoErro(tempo, solucao_aproximada, solucao_exata, metodo, ordem):
    
    #print("len tempo: ", len(tempo))
    #print("len solucao_aproximada: ", len(solucao_aproximada))
    #print("len solucao_exata: ", len(solucao_exata))
    
    arr_erro = np.zeros((len(tempo),), dtype=np.float128)
    
    #gera array com a comparacao do erro
    for k in range (0, num_particoes):
        arr_erro[k] = abs(solucao_aproximada[k] - solucao_exata[k])
        
    
     #plota grafico da função
    plt.plot(
        tempo, arr_erro, 'b--'
        )
    #plt.ylabel(u"Erro f(x) - P(x)") #esse 'u' antes da string é pra converter o texto pra unicode
    #plt.xlabel(u"Valor de x, " + str(len(tempo)) + u" partições")
    plt.xlabel(str(len(tempo)) + u" partições, N = " + str(ordem))
    
    #legendas do grafico
    se_line = mlines.Line2D([], [], color='blue', marker='', markersize=0, label=u'Erro')
    
    plt.legend(handles=[se_line])
    
    plt.title(u"Erro relativo do metodo "+metodo)
    
    #plt.axis([0, 50, 0, 100])
    plt.show() 
    


#configura grafico
num_particoes = 100
N = 8
a = -1.0
b = 1.0


#gera os pontos pra interpolacao
pontos = gerarPontos_a(a, b, N)


#inicia os arrays do numpy com 128bits de precisao
particoes = np.zeros((num_particoes,), dtype=np.float128)
sol_exata = np.zeros((num_particoes,), dtype=np.float128)
sol_aprox = np.zeros((num_particoes,), dtype=np.float128)



#gera array com a solucao exata e particoes do grafico
dp = (b-a) / num_particoes
for k in range (0, num_particoes):
    x = a + dp * k
    particoes[k] = x
    sol_exata[k] = F(x)
    
    

#gera array com pontos pra interpolacao
num_pontos = len(pontos)
y = np.zeros((num_pontos,), dtype=np.float128)
for k in range (0, num_pontos):
    y[k] = F(pontos[k])  
    
'''
print(particoes)
print(pontos)
print(y)
'''

'''
print("")
print("------")
print("")

#inicio = Utils().getTime()
#Newton().interpolacao(pontos, vy);
coeficientes = Newton().interpolacaoSlide(num_pontos, pontos, y);

#gera array com a solucao aproximada
for k in range (0, num_particoes):
    #sol_aprox[k] = Newton().avaliarPonto(particoes[k])
    sol_aprox[k] = Newton().avaliarPontoSlide(num_pontos, particoes[k], pontos, coeficientes)
    
#fim  = Utils().getTime()
#Utils().imprimeDiferencaTempo(inicio, fim)
erroNorma = Utils().distanciaMaximo(sol_exata, sol_aprox)
print("Erro (Norma do Maximo): " + repr(erroNorma))


gerarGrafico(particoes, sol_aprox, sol_exata, "Newton", N)
gerarGraficoErro(particoes, sol_aprox, sol_exata, "Newton", N)
'''

print("")
print("------")
print("")

#inicio = Utils().getTime()

#gera array com a solucao aproximada
for k in range (0, num_particoes):
    sol_aprox[k] = Lagrange().interpolacao(num_pontos, pontos, y, particoes[k]);
    
#fim  = Utils().getTime()
#Utils().imprimeDiferencaTempo(inicio, fim)
erroNorma = Utils().distanciaMaximo(sol_exata, sol_aprox)
print("Erro (Norma do Maximo): " + repr(erroNorma))

gerarGrafico(particoes, sol_aprox, sol_exata, "Lagrange", N)
gerarGraficoErro(particoes, sol_aprox, sol_exata, "Lagrange", N)



print("")
print("------")
print("")

#inicio = Utils().getTime()

#gera array com a solucao aproximada
for k in range (0, num_particoes):
    sol_aprox[k] = Linear().interpolacao(num_pontos, pontos, y, particoes[k]);
    
#fim  = Utils().getTime()
#Utils().imprimeDiferencaTempo(inicio, fim)
erroNorma = Utils().distanciaMaximo(sol_exata, sol_aprox)
print("Erro (Norma do Maximo): " + repr(erroNorma))

gerarGrafico(particoes, sol_aprox, sol_exata, "Linear", 1)
gerarGraficoErro(particoes, sol_aprox, sol_exata, "Linear", 1)
