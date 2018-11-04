# -*- coding: utf-8 -*-
"""
Created on Fri Nov  2 14:12:16 2018

@author:  Thiago Almeida
"""

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.pyplot import figure
import numpy as np

from MinimosQuadrados import MinimosQuadrados



#define o tamanho dos graficos
figure(num=None, figsize=(8, 6), dpi=72, facecolor='w', edgecolor='k')


###### Leitura do arquivo de dados #######

arquivo = open("dados.txt", 'r');

x = np.empty((0,), dtype=np.float128)
y = np.empty((0,), dtype=np.float128)


for linha in arquivo:
    #quebra a linha em uma lista, separando pela tabulacao
    dados = linha.split('\t')
    
    #tratamento das strings pra conversao pra numero
    a = dados[0].replace(",", ".")
    b = dados[1].replace("\n", "").replace(",", ".")
    
    #adiciona ao array
    x = np.append(x, np.float128(a))
    y = np.append(y, np.float128(b))


print("len(x)", len(x))
print("len(y)", len(y))


def gerarGrafico(x, y, x_aprox, y_aprox, ordem):
    
    #plota pontos
    for k in range (0, len(x)):
        plt.plot(x[k], y[k], 'yo')
            
    #plota grafico da função
    plt.plot(
        x_aprox, y_aprox, 'b-'  
        )
    plt.ylabel(u"Cotação do dia em reais") #esse 'u' antes da string é pra converter o texto pra unicode
    plt.xlabel(u"dias dos ultimos 2 anos")
    
    #legendas do grafico
    se_line = mlines.Line2D([], [], color='blue', marker='', markersize=0, label=u'Aproximação (N='+str(ordem)+')')
    ac_line = mlines.Line2D([], [], color='yellow', marker='', markersize=0, label=u'Dados')
    
    plt.legend(handles=[se_line, ac_line])
    
    plt.title(u"Cotações da Petrobrás nos ultimos 2 anos", )
    
    #plt.axis([0, 50, 0, 100])
    plt.show()
    
###########   
    
N = 8
num_pontos = 100

###########

#pontos do slide, que geram o grafico da pagina 31
#x = np.array([-1, 0, 1, 2], dtype=np.float128)
#y = np.array([0, -1, 0, 7], dtype=np.float128)


#calcula minimos quadrados
coeficientes = MinimosQuadrados().executar(x, y, N)

#cria vetores pra interpolacao
z = np.zeros((num_pontos+1,), dtype=np.float128)
pontos = np.zeros((num_pontos,), dtype=np.float128)
particoes = np.zeros((num_pontos,), dtype=np.float128)

#gera array com a solucao exata e particoes do grafico
a = min(x)
b = max(x)
dp = (b-a) / num_pontos
particoes = np.arange(a, b, dp)
particoes = np.append(particoes, b)

#calcula pontos interpolados
for k in range (0, num_pontos+1):
    z[k] = MinimosQuadrados().interpolaCoeficientes(coeficientes, N, particoes[k])
    
    
gerarGrafico(x, y, particoes, z, N)