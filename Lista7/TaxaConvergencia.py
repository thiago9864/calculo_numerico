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

figure(num=None, figsize=(8, 6), dpi=72, facecolor='w', edgecolor='k')

#numero de elementos
x = [4, 16, 64, 256, 1024]

#intervalo
a = 0.0
b = 1.0

####### Dados dos erros ########



#montagem do vetor: [erro com 4 elementos, erro com 5 elementos, ... , erro com 1024 elementos]
n1 = [9, 8, 7, 6, 5]
n2 = [8, 7, 6, 5, 4]
n3 = [7, 6, 5, 4, 3]
n4 = [6, 5, 4, 3, 2]
n5 = [5, 4, 3, 2, 1]



####### tira log dos dados ########

#faz -log(h)
for i in range(5):
    h = (b-a)/x[i]
    x[i] = -np.log(h)

#faz log(erro)    
for i in range(5):
    n1[i] = np.log(n1[i])
    n2[i] = np.log(n2[i])
    n3[i] = np.log(n3[i])
    n4[i] = np.log(n4[i])
    n5[i] = np.log(n5[i])

####### Criacao do grafico ########

print(n5)

#plota grafico do N=1
plt.plot(
    x, n1, 'r-'
    )
plt.plot(
    x, n2, 'g-'
    )
plt.plot(
    x, n3, 'b-'
    )
plt.plot(
    x, n4, 'm-'
    )
plt.plot(
    x, n5, 'k-'
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