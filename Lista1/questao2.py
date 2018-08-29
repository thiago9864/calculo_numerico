# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import math as m
from array import array


nParticoes = 4**2

#Solicitar o valor de E
#testes com E = 0.1, 0.01, 0.001, 0.0001
E = 0.001


#Coeficientes
a = 0
b = 1
Nel = nParticoes - 1
h = (b-a)/float(Nel)
c2 = (m.exp((-1)/m.sqrt(E)) - 1) / (m.exp(1/m.sqrt(E)) - m.exp((-1)/m.sqrt(E))) 
c1 = - 1 - c2

#Definicao da funcao
def F(x): 
    return c1 * m.exp(- x / m.sqrt (E)) +  c2 * m.exp(x / m.sqrt (E)) + 1

def discretizacao(a, b, c):
    return (a * u_ant + b * u + u_prox * c == h**2)



#monta matriz
ordem = nParticoes - 2
ordem_is = ordem - 1
ds = []
dp = []
di = []
d = []
arr_que_sei_la = []  
arr_tempo = []
arr_solucao_exata = []

#preenche diagonal superior e inferior
for i in range(ordem_is):
    ds.append(E * -1)
    di.append(E * -1)
    
#preenche diagonal principal
for i in range(ordem):
    dp.append(2 * E * + h**2)
    d.append(-h**2)
    dh = h * i
    arr_tempo.append(dh)
    arr_solucao_exata.append(F(dh))
    

#Algoritmo de Thomas
def Thomas(a, b, c, d):
    ''' Resolve Ax = d onde A e uma matriz tridiagonal composta pelos vetores a, b, c
		a - subdiagonal
        	b - diagonal principal
		c - superdiagonal.
	Retorna x
    '''
    #print('diagonais: ', a, b, c, d)
    n = len(d) # len(d) == len(b)
    c_ = [ c[0] / b[0] ]
    d_ = [ d[0] / b[0] ]
    
    for i in range(1, n):
        aux = b[i] - c_[i-1]*a[i-1]
        if i < n-1:
            c_.append( c[i] / aux )
        d_.append( (d[i] - d_[i-1]*a[i-1])/aux )
    
    # Substituicao de volta
    x = [d_[-1]]
    for i in range(n-2, -1, -1):
        x = [ d_[i] - c_[i]*x[0] ] + x
    
    return x
    
  
arr_que_sei_la = Thomas(di, dp, ds, d)

#plota grafico
plt.plot(
    arr_tempo, arr_solucao_exata, 'b--',
    arr_tempo, arr_que_sei_la, 'r-'
    )
plt.ylabel(u"Temperatura (ºC)") #esse 'u' antes da string é pra converter o texto pra unicode
plt.xlabel(u"Tempo (seg), " + str(nParticoes) + u" partições")

#legendas do grafico
se_line = mlines.Line2D([], [], color='blue', marker='', markersize=0, label=u'Solução Exata')
ac_line = mlines.Line2D([], [], color='red', marker='', markersize=0, label=u'Alguma coisa')

plt.legend(handles=[se_line, ac_line])

plt.title("Metodos de Resolucao")

#plt.axis([0, 50, 0, 100])
plt.show()