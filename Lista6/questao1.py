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
from MetodosIntegracao import QuadraturaGauss
from MinimosQuadrados import MinimosQuadrados


#define o tamanho dos graficos
figure(num=None, figsize=(8, 6), dpi=72, facecolor='w', edgecolor='k')


###### Testes #######
'''
#funcao do slide
def f(x):
    return x**2 + 2.0*x + 2    
    r#eturn np.exp(x) * np.cos(x)

#intervalo
a = 0.0
b = 1.2
num_divisoes = 6 * 100
num_divisoes2 = num_divisoes+1
num_divisoes3 = num_divisoes+1

num_pontos_gauss = 2

print("num_divisoes", num_divisoes)
print("num_pontos_gauss", num_pontos_gauss)

#valor_exato = 14.0 / 3.0 #4.666666666666667 (a, b) = [-1, 1]
valor_exato = 4.416 #(a, b) = [0, 1.2]

#valor_exato = 1.648774427347253


repetido_retangulo = Repetidos().Retangulo(a, b, num_divisoes, f)
repetido_pontomedio = Repetidos().PontoMedio(a, b, num_divisoes, f)
repetido_trapezio = Repetidos().Trapezio(a, b, num_divisoes, f)
repetido_simpsom13 = Repetidos().Simpsom13(a, b, num_divisoes2, f)
repetido_simpsom38 = Repetidos().Simpsom38(a, b, num_divisoes3, f)
repetido_gauss = QuadraturaGauss().Gauss(a, b, num_pontos_gauss, f)

print("Valor exato", valor_exato)
print("NewtonCotes().Retangulo", NewtonCotes().Retangulo(a, b, f))
print("NewtonCotes().PontoMedio", NewtonCotes().PontoMedio(a, b, f))
print("NewtonCotes().Trapezio", NewtonCotes().Trapezio(a, b, f))
print("NewtonCotes().Simpsom13", NewtonCotes().Simpsom13(a, b, f))
print("NewtonCotes().Simpsom38", NewtonCotes().Simpsom38(a, b, f))
print("--- repetidos ---")
print("Repetidos().Retangulo", repetido_retangulo, abs(valor_exato - repetido_retangulo))
print("Repetidos().PontoMedio", repetido_pontomedio, abs(valor_exato - repetido_pontomedio))
print("Repetidos().Trapezio", repetido_trapezio, abs(valor_exato - repetido_trapezio))
print("Repetidos().Simpsom13", repetido_simpsom13, abs(valor_exato - repetido_simpsom13))
print("Repetidos().Simpsom38", repetido_simpsom38, abs(valor_exato - repetido_simpsom38))
print("--- Gauss ---")
print("QuadraturaGauss().Gauss", repetido_gauss, abs(valor_exato - repetido_gauss))
'''

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

def f(x):
    global c1, c2, E
    return u(c1, c2, E, x)
    #return x**4 - 5.0*x


#intervalo de integracao
a = 0.0
b = 1.0

#pontos do Gauss
N = 10

#pontos pra gerar o grafico
num_pontos = 100


###### GRAFICOS ######


def gerarGrafico(x, y, x_aprox, y_aprox, ordem):
    
    print(len(x), len(y))
    print(len(x_aprox), len(y_aprox))
            
    #plota grafico da função
    plt.plot(
        x, y, 'b-',
        x_aprox, y_aprox, 'r-' 
        )
    plt.ylabel(u"f(x)") #esse 'u' antes da string é pra converter o texto pra unicode
    plt.xlabel(u"x")
    
    #legendas do grafico
    se_line = mlines.Line2D([], [], color='blue', marker='', markersize=0, label=u'Função Exata')
    ac_line = mlines.Line2D([], [], color='red', marker='', markersize=0, label=u'Aproximação (N='+str(ordem)+')')
    
    plt.legend(handles=[se_line, ac_line], loc='lower center')
    
    '''Posicoes da legenda 
        upper right
        upper left
        lower left
        lower right
        right
        center left
        center right
        lower center
        upper center
        center
    '''
    
    plt.title(u"Função exata e interpolação", )
    
    plt.show()
    
    

###### LETRA A ######


#cria vetores pra interpolacao
y = np.zeros((num_pontos,), dtype=np.float128)
z = np.zeros((num_pontos,), dtype=np.float128)
pontos = np.zeros((num_pontos,), dtype=np.float128)
particoes = np.zeros((num_pontos,), dtype=np.float128)


#gera array com a solucao exata e particoes do grafico
dp = (b-a) / (num_pontos-1)
soma = a
for k in range (0, num_pontos):
    particoes[k] = soma
    soma += dp
    
#calcula solucao exata
for k in range (0, num_pontos):
    y[k] = f(particoes[k])

#calcula minimos quadrados
coeficientes = MinimosQuadrados().executar(a, b, N, f)

#calcula pontos interpolados
for k in range (0, num_pontos):
    z[k] = MinimosQuadrados().interpolaCoeficientes(coeficientes, N, particoes[k])

   
gerarGrafico(particoes, y, particoes, z, N)

#print("coeficientes", coeficientes)
#print("particoes", particoes)
#print("y", y)
#print("z", z)