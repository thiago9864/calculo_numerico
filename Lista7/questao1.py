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

from Erros import Erro
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
    
# Funcao u(x)    
E = 0.001
t = 1.0 / np.sqrt(E)
c1 = -1.0
c2 = (np.exp(-t) - 1.0) / (np.exp(t) - np.exp(-t))

def f(x):
    global c1, c2, E
    return u(c1, c2, E, x)
    #return x**4 - 5.0*x


###### Parametros #######

#intervalo
a = 0
b = 1.0

#contorno
y1 = 0
y2 = 0

#configuracoes
elementos = 4
grau_polinomio = 2
form_fraca = 1


###### Inicia Elementos Finitos #######

elementos_finitos = ElementosFinitos(a, b, y1, y2, elementos, grau_polinomio, form_fraca, dataType)


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

q = (b-a) / (particoes_exata-1)
d = a

for k in range(particoes_exata):
    x[k] = d
    d += q
            
for k in range(particoes_exata):
    y[k] = f_slide(x[k])
  
print("x")
print(x)
print("y")
print(y)


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
    se_line = mlines.Line2D([], [], color='blue', marker='', markersize=0, label=u'Função Exata')
    ac_line = mlines.Line2D([], [], color='red', marker='', markersize=0, label=u'Aproximação ('+str(elementos)+' elem, pol ordem '+str(ordem)+')')
    
    #plt.legend(handles=[se_line, ac_line], loc='lower center')
    
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
    #plt.axis([0, 1, 0, 1.05])
    
    plt.show()

gerarGrafico(x, y, xu, U, elementos, grau_polinomio)