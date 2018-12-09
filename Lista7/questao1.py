# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 11:12:42 2018

@author: thiagoalmeida
"""

#import platform
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.pyplot import figure
import numpy as np

from ElementosFinitos import ElementosFinitos

'''
muda pra float128 se for windows
'''
dataType = np.float64
#if(platform.system() != 'Windows'):
#    dtype = np.float128
    
try:
    dataType = np.float128
    print("Sistema tem suporte a float128")
except:
    dataType = np.float64
    print("Sistema não tem suporte a float128, usando float64")


#define o tamanho dos graficos
figure(num=None, figsize=(8, 6), dpi=72, facecolor='w', edgecolor='k')

###### Funcao de teste da apostila #######

def f_slide(x):
    return (-(x**2) / 2.0) + (x/2.0)
    
###### Funcoes da lista #######


def u(c1, c2, E, x):
    t = x / np.sqrt(E)
    return c1 * np.exp(-t) + c2 * np.exp(t) + 1.0
    
    
###### Definicoes #######
       
E = 10**(-4)

t = 1.0 / np.sqrt(E)
c1 = -1.0
c2 = (np.exp(-t) - 1.0) / (np.exp(t) - np.exp(-t))

def f(x):
    global c1, c2, E
    return u(c1, c2, E, x)
    #return x**4 - 5.0*x

print("E:", E)

###### Parametros #######

#intervalo
a = 0
b = 1.0

#contorno
y1 = 0
y2 = 0

#configuracoes
elementos = 42
grau_polinomio = 2


###### Inicia Elementos Finitos #######

elementos_finitos = ElementosFinitos(a, b, y1, y2, elementos, grau_polinomio, E, dataType)


K = elementos_finitos.matrizLocal()
F = elementos_finitos.vetorLocal()

print("K")
print(K)
print("\nF")
print(F)

Kr = elementos_finitos.matrizRigidez()
Fr = elementos_finitos.vetorForca()

print("\nKr")
print(Kr)
print("\nFr")
print(Fr)

U = elementos_finitos.calcular()

print("\nU")
print(U)

len_num = len(U)
xu = np.zeros((len_num,), dtype=dataType)
hu = (b-a) / (len_num-1)
d = a
for k in range(len_num):
    xu[k] = d
    d += hu

###### Cria dados da funcao exata #######


particoes_exata = 100


x = np.zeros((particoes_exata,), dtype=dataType)
y = np.zeros((particoes_exata,), dtype=dataType)
yu = np.zeros((len_num,), dtype=dataType)

q = (b-a) / (particoes_exata-1)
d = a

#cria vetor x de particoes pra solucao exata
for k in range(particoes_exata):
    x[k] = d
    d += q

#cria vetor y da funcao exata nos pontos x do vetor criado com 100 particoes        
for k in range(particoes_exata):
    y[k] = f(x[k])#lista
    #y[k] = f_slide(x[k])#slide
    
#cria vetor y da funcao exata nos pontos x do vetor dos elementos finitos       
for k in range(len_num):
    yu[k] = f(xu[k])#lista
    #yu[k] = f_slide(xu[k])#slide
    
'''  
print("x")
print(x)
print("y")
print(y)
'''

###### GRAFICOS ######


def gerarGrafico(x, y, x_aprox, y_aprox, elementos, ordem):
    
    #print(len(x), len(y))
    #print(len(x_aprox), len(y_aprox))
            
    #plota grafico da função
    plt.plot(
        x, y, 'b-'
        )
        
    plt.plot(
        x_aprox, y_aprox, 'r-' 
        )
    plt.ylabel(u"f(x)") #esse 'u' antes da string é pra converter o texto pra unicode
    plt.xlabel(u"x")
    
    
    #legendas do grafico    
    se_line = mlines.Line2D([], [], color='blue', marker='', markersize=0, label=u'Solução Exata')
    ac_line = mlines.Line2D([], [], color='red', marker='', markersize=0, label=u'Aproximação ('+str(elementos)+' elementos, n = '+str(ordem)+')')
    
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
    
    plt.title(u"Solução Exata e Aproximação", )
    
    #muda os limites dos eixos
    #por: x_inicial, x_final, y_inicial, y_final
    #plt.axis([0, 1, 0, 1.05])
    
    plt.show()

gerarGrafico(x, y, xu, U, elementos, grau_polinomio)