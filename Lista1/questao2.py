# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import math as m
from array import array

#Definicao da funcao
def F(x, c1, c2, E): 
    return c1 * m.exp(- x / m.sqrt (E)) +  c2 * m.exp(x / m.sqrt (E)) + 1.0

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
    

def geraGraficoFuncao(epsilon, indice, plotar):
    
    nParticoes = 4**indice
    
    nParticoes = 53

    #Solicitar o valor de E
    #testes com E = 0.1, 0.01, 0.001, 0.0001
    E = epsilon
    
    #Coeficientes
    a = 0
    b = 1.0
    Nel = nParticoes
    h = (b-a)/float(Nel)
    c2 = (m.exp((-1)/m.sqrt(E)) - 1) / (m.exp(1/m.sqrt(E)) - m.exp((-1)/m.sqrt(E))) 
    c1 = - 1 - c2
    
    h2 = h**2.0
    cdp = 2 * E + h2

    #monta matriz
    ordem = nParticoes - 1
    ordem_is = ordem - 1
    ds = []
    dp = []
    di = []
    d = []
    arr_solucao_aproximada = []  
    arr_tempo_exata = []
    arr_tempo = []
    arr_solucao_exata = []
    arr_erro = []
    
    #preenche diagonal superior e inferior
    for i in range(ordem_is):
        ds.append(E * -1.0)
        di.append(E * -1.0)
        
    #preenche diagonal principal
    for i in range(ordem):
        dp.append(cdp)
        d.append(h2)
        
    #calcula tempo e solucao exata
    particoes_exata = 100
    he = (b-a)/float(particoes_exata-1)
    for i in range(particoes_exata):
        dhe = he * i
        arr_tempo_exata.append(dhe)
        arr_solucao_exata.append(F(dhe, c1, c2, E))  
        
    #calcula tempo da solucao aproximada
    for i in range(nParticoes+1):
        dh = h * i
        arr_tempo.append(dh)
        
    #resolve sistema para solucao aproximada
    arr_solucao_aproximada = Thomas(di, dp, ds, d)
    
    #insere intervalos de contorno
    arr_solucao_aproximada.insert(0, 0)
    arr_solucao_aproximada.append(0)
    
    #calcula erro
    for i in range(nParticoes+1):
        err = m.fabs(arr_solucao_exata[i] - arr_solucao_aproximada[i]);
        arr_erro.append(err)
    
    if(plotar):
        print("plotar")
        #plota grafico da função
        plt.plot(
            arr_tempo_exata, arr_solucao_exata, 'b--' 
            )
        
        plt.plot(
            arr_tempo, arr_solucao_aproximada, 'r-'    
            )
        
        plt.ylabel(u"Valor de u(h)") #esse 'u' antes da string é pra converter o texto pra unicode
        plt.xlabel(u"Valor de h, " + str(nParticoes) + u" partições")
        
        #legendas do grafico
        se_line = mlines.Line2D([], [], color='blue', marker='', markersize=0, label=u'Solução Exata')
        ac_line = mlines.Line2D([], [], color='red', marker='', markersize=0, label=u'Solução Aprox.')
        
        plt.legend(handles=[se_line, ac_line])
        
        plt.title("Metodos de Resolucao")
        
        #plt.axis([0, 50, 0, 100])
        plt.show()    
    
    return [arr_tempo, arr_erro]
    
#executa o codigo
    
geraGraficoFuncao(0.001, 2, 1)
    
'''
lst_epsilon = [0.1, 0.01, 0.001, 0.0001]
lst_indices = [1, 2, 3, 4, 5]

lst_erro = geraGraficoFuncao(0.001, 2, 0)

plt.plot(
        lst_erro[0], lst_erro[1], 'g-'  
        )
        
er_line = mlines.Line2D([], [], color='green', marker='', markersize=0, label=u'Erro.')

plt.legend(handles=[er_line])
plt.title("Erro")

plt.show()
'''