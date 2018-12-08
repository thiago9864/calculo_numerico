#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  8 17:53:04 2018

@author: thiagoalmeida
"""

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.pyplot import figure
import numpy as np

from Erros import Erro
from ElementosFinitos import ElementosFinitos

figure(num=None, figsize=(8, 6), dpi=72, facecolor='w', edgecolor='k')

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


#intervalo
a = 0.0
b = 1.0

#contorno
y1 = 0.0
y2 = 0.0

###### Definicoes #######
       
E = 10**(-3)

t = 1.0 / np.sqrt(E)
c1 = -1.0
c2 = (np.exp(-t) - 1.0) / (np.exp(t) - np.exp(-t))

def u(c1, c2, E, x):
    t = x / np.sqrt(E)
    return c1 * np.exp(-t) + c2 * np.exp(t) + 1.0

def f(x):
    global c1, c2, E
    return u(c1, c2, E, x)
    #return x**4 - 5.0*x


####### Dados dos erros ########

x = np.zeros((5,), dtype=dataType)
n = np.zeros((5,5), dtype=dataType)

for i in range(1, 4):#elementos
    for j in range(1, 6):#grau do polinomio
               
        elementos = int(4.0**i)
        grau_polinomio = int(j)
        
        print("N: "+str(grau_polinomio)+", Elementos: " + str(elementos))
        
        elementos_finitos = ElementosFinitos(a, b, y1, y2, elementos, grau_polinomio, E, dataType)
        U = elementos_finitos.calcular()
        
        len_num = len(U)
        
        xu = np.zeros((len_num,), dtype=dataType)
        yu = np.zeros((len_num,), dtype=dataType)
        
        hu = (b-a) / (len_num-1)
        d = a
        for k in range(len_num):
            xu[k] = d
            d += hu
    
        #cria vetor y da funcao exata nos pontos x do vetor dos elementos finitos       
        for k in range(len_num):
            yu[k] = f(xu[k])#lista
            #yu[k] = f_slide(xu[k])#slide
        
        #print(xu)
        #print(U)
        #print(yu)
        print("-------")
        

        ########### Calcular o erro aqui ###########
        
        erro = 0
        
        ############################################
        
        n[j-1][i-1] = np.log(erro)
        
    x[i-1] = (b-a) / int(4.0**i)

'''
#montagem do vetor: [erro com 4 elementos, erro com 5 elementos, ... , erro com 1024 elementos]
n1 = [9, 8, 7, 6, 5]
n2 = [8, 7, 6, 5, 4]
n3 = [7, 6, 5, 4, 3]
n4 = [6, 5, 4, 3, 2]
n5 = [5, 4, 3, 2, 1]
'''
print(x)
print(n)

####### tira log dos dados ########

#faz -log(h)
for i in range(5):
    x[i] = -np.log(x[i])

#faz log(erro) 
'''
for i in range(5):
    n1[i] = np.log(n1[i])
    n2[i] = np.log(n2[i])
    n3[i] = np.log(n3[i])
    n4[i] = np.log(n4[i])
    n5[i] = np.log(n5[i])
'''


####### Criacao do grafico ########


#plota grafico do N=1
plt.plot(
    x, n[0], 'r-'
    )
plt.plot(
    x, n[1], 'g-'
    )
plt.plot(
    x, n[2], 'b-'
    )
plt.plot(
    x, n[3], 'm-'
    )
plt.plot(
    x, n[4], 'k-'
    )

'''Cores dos graficos
b: blue
g: green
r: red
c: cyan
m: magenta
y: yellow
k: black
w: white
'''

leg_n1 = mlines.Line2D([], [], color='red', marker='', markersize=0, label=u'N=1')
leg_n2 = mlines.Line2D([], [], color='green', marker='', markersize=0, label=u'N=2')
leg_n3 = mlines.Line2D([], [], color='blue', marker='', markersize=0, label=u'N=3')
leg_n4 = mlines.Line2D([], [], color='magenta', marker='', markersize=0, label=u'N=4')
leg_n5 = mlines.Line2D([], [], color='black', marker='', markersize=0, label=u'N=5')

plt.legend(handles=[leg_n1, leg_n2, leg_n3, leg_n4, leg_n5], loc='upper right')
    
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



plt.ylabel(u"log(erro)") #esse 'u' antes da string é pra converter o texto pra unicode
plt.xlabel(u"-log(h)")
plt.title(u"Taxa de Convergência", )

#muda os limites dos eixos
#por: x_inicial, x_final, y_inicial, y_final
#plt.axis([-5, 5, -5, 5])

plt.show()