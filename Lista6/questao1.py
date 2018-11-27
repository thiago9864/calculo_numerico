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
from Erros import Erro


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
N = 5

#pontos pra gerar o grafico
num_pontos = 100


###### GRAFICOS ######


def gerarGrafico(x, y, x_aprox, y_aprox, ordem):
    
    #print(len(x), len(y))
    #print(len(x_aprox), len(y_aprox))
            
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
    
    #muda os limites dos eixos
    #por: x_inicial, x_final, y_inicial, y_final
    plt.axis([0, 1, 0, 1.05])
    
    plt.show()
    

print("\n###### LETRA A ######\n")


#cria vetores pra interpolacao
y = np.zeros((num_pontos,), dtype=np.float64)
z = np.zeros((num_pontos,), dtype=np.float64)
pontos = np.zeros((num_pontos,), dtype=np.float64)
particoes = np.zeros((num_pontos,), dtype=np.float64)


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
coeficientes_a = MinimosQuadrados().executar(a, b, N, f)

#calcula pontos interpolados
for k in range (0, num_pontos):
    z[k] = MinimosQuadrados().interpolaCoeficientes(coeficientes_a, N, particoes[k])

   
gerarGrafico(particoes, y, particoes, z, N)

#print("coeficientes", coeficientes)
#print("particoes", particoes)
#print("y", y)
#print("z", z)


print("\n###### LETRA B ######\n")

# n nessa letra fica fixado em 1
Nb = 1


#gera os 25 intervalos
num_subintervalos = 25
subintervalos = [];
dp_s = (b-a) / num_subintervalos

ini = 0.0
fim = 0
pi = 0
pf = 0
for k in range (0, num_subintervalos):
    #fim dos indices de pontos
    fim = ini + dp_s
    #fim da divisao de intervalos
    pf = pi + 4
    #monta vetor de indices do vetor de particoes
    part = np.arange(pi, pf, 1)
    #monta entrada de intervalo na lista de intervalos
    subintervalos.append([ini, fim, list(part)])
    #atualiza pro proximo loop
    ini = fim
    pi = pf


#print(len(subintervalos), subintervalos)
erro_curva_b = 0

#calcula os pontos pra cada um dos subintervalos
for k in range (0, num_subintervalos):
    _a = subintervalos[k][0]
    _b = subintervalos[k][1]
    rg = subintervalos[k][2]
    
    #calcula minimos quadrados
    coeficientes_b = MinimosQuadrados().executar(_a, _b, Nb, f)
    

    #calcula pontos interpolados
    for k in rg:
        z[k] = MinimosQuadrados().interpolaCoeficientes(coeficientes_b, Nb, particoes[k])
        
    #indices pra cortar a lista de  particoes e mandar as 4 necessarias    
    p = rg[0]
    q = rg[3]+1
    
    #calcula erro de cada subdivisao
    erro_curva_b += Erro().erroDaNorma(4,_a,_b,f, coeficientes_b, Nb, particoes[p:q])

gerarGrafico(particoes, y, particoes, z, Nb)


print("\n###### LETRA C ######\n")


erro_curva_a = Erro().erroDaNorma(num_pontos,a,b,f, coeficientes_a, N, particoes)
print ("Erro da letra A, com N=" + str(N) + ": " + repr(erro_curva_a))


print ("Erro da letra B, com N=" + str(Nb) + ": " + repr(erro_curva_b))