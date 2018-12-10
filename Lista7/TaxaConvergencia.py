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

for i in range(1, 5):#elementos
    for j in range(1, 6):#grau do polinomio
               
        elementos = int(4.0**i)
        grau_polinomio = int(j)
        
        print("N: "+str(grau_polinomio)+", Elementos: " + str(elementos))
        
        elementos_finitos = ElementosFinitos(a, b, y1, y2, elementos, grau_polinomio, E, dataType)
        U = elementos_finitos.calcular()
        
        len_num = len(U)
        
        #calcula funcao exata
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
        
        erro = Erro().calculaErro(yu, U)        
        
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

####### tira log dos dados ########

for i in range(5):
    x[i] = -np.log(x[i])

print(x)
print(n)


#Matriz de dados prontos
'''   
x = [0.25,     0.0625,   0.015625, 0,       0,      ]
n = [[-1.43493175e+00, -2.77528630e+00, -5.56543749e+00, -8.35931975e+00,   0.00000000e+00]
     [-2.45211323e+00, -5.64960746e+00, -1.09238778e+01, -1.60678043e+01,   0.00000000e+00]
     [-2.37346130e+00, -1.10547675e+00, -1.36817671e-01,  1.40952962e-01,   0.00000000e+00]
     [-1.94081421e+00, -6.60388993e-01, -7.53548626e-02,  8.30020979e-02,   0.00000000e+00]
     [-1.68315286e+00, -4.64059450e-01, -1.56213424e-02,  1.09499787e-01,   0.00000000e+00]]
'''
erro_dif_finitas = [0.01514269848800638, 0.03590183209403578, 0.003660693293331358, 0.0002335478042604083, 1.4615914875454955e-05]

for i in range(5):
    erro_dif_finitas[i] = np.log(erro_dif_finitas[i])
    
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
plt.plot(
    x, erro_dif_finitas, 'c-'
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
leg_nd = mlines.Line2D([], [], color='cyan', marker='', markersize=0, label=u'Dif. Finitas')

plt.legend(handles=[leg_n1, leg_n2, leg_n3, leg_n4, leg_n5, leg_nd], loc='lower left')
    
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