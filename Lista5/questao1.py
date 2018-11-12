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

#tive que fazer isso pra rodar no Linux
conteudo = arquivo.readlines();
string = ""
for i in conteudo:
   string += i


#prestar atenção nisso, se for windows deve ter q mudar esse if
if(string.find('\r') >= 0):
    #cai nesse se for Linux
    linhas = string.split('\r')
else:
    #cai nesse se for mac
    linhas = string.split('\n')
    
#print(linhas)
    

for linha in linhas:
    #quebra a linha em uma lista, separando pela tabulacao
    
    dados = linha.split('\t')
    
    #tratamento das strings pra conversao pra numero
    a = dados[0].replace(',', '.')
    b = dados[1].replace("\n", "").replace("\r", "").replace(',', '.')
    
    #adiciona ao array
    x = np.append(x, np.float128(a))
    y = np.append(y, np.float128(b))


print("len(x)", len(x))
print("len(y)", len(y))


def gerarGrafico(x, y, x_aprox, y_aprox, ordem):
    
    print(len(x_aprox), len(y_aprox))
    
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
    
    plt.legend(handles=[se_line, ac_line], loc='upper left')
    
    plt.title(u"Cotações da Petrobrás nos ultimos 2 anos", )
    
    plt.show()
    
    
def gerarGraficoPrevisao(x, y, x_aprox, y_aprox, x_prev, y_prev, ordem):
    
    #plota pontos
    for k in range (0, len(x)):
        plt.plot(x[k], y[k], 'yo')
            
    #plota grafico da função
    plt.plot(
        x_aprox, y_aprox, 'b-'
        )
        
    print(len(x_prev), len(y_prev))
        
    plt.plot(
        x_prev, y_prev, 'r-' 
        )
        
    plt.ylabel(u"Cotação do dia em reais") #esse 'u' antes da string é pra converter o texto pra unicode
    plt.xlabel(u"dias dos ultimos 2 anos + 100 dias futuros")
    
    #legendas do grafico
    se_line = mlines.Line2D([], [], color='blue', marker='', markersize=0, label=u'Aproximação (N='+str(ordem)+')')
    ac_line = mlines.Line2D([], [], color='yellow', marker='', markersize=0, label=u'Dados')
    pr_line = mlines.Line2D([], [], color='red', marker='', markersize=0, label=u'Previsão')
    
    plt.legend(handles=[se_line, ac_line, pr_line], loc='upper left')
    
    plt.title(u"Previsão das cotações nos próximos 100 dias", )
    
    #plt.axis([0, 700, 10, 30])
    plt.show()
    
###########   
    
N = 95
num_pontos = 100

###########

#pontos do slide, que geram o grafico da pagina 31
#x = np.array([-1, 0, 1, 2], dtype=np.float128)
#y = np.array([0, -1, 0, 7], dtype=np.float128)


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


#calcula minimos quadrados
coeficientes = MinimosQuadrados().executar(x, y, N)

#calcula pontos interpolados
for k in range (0, num_pontos+1):
    z[k] = MinimosQuadrados().interpolaCoeficientes(coeficientes, N, particoes[k])

#calcula o r    
r = MinimosQuadrados().calculaResiduo(y, x, N, coeficientes)
   
gerarGrafico(x, y, particoes, z, N)

print("coeficientes", coeficientes)
print("particoes", particoes)
print("z", z)
print("N = "+repr(N)+", r = " + repr(r));

'''
print("")
print("######## Calcula melhor possivel ##########")
print("")

##### calcula todos os r para cada n de 0 a 100 #####

res = np.zeros((num_pontos,), dtype=np.float128)

#print(res)

for ni in range (1, num_pontos):
    
    c2 = MinimosQuadrados().executar(x, y, ni)
    
    #calcula pontos interpolados
    for k in range (0, num_pontos+1):
        z[k] = MinimosQuadrados().interpolaCoeficientes(c2, ni, particoes[k])
        
    #calcula o r
    r = MinimosQuadrados().calculaResiduo(y, x, ni, c2)
    res[ni] = r
    
print("res (valores de r para cada ordem)", res)

#pega o maior valor
maior = res[0]
melhor_indice = 0
for k in range (0, len(res)):
    if(res[k] > maior):
        maior = res[k]
        melhor_indice = k

print("")        
print("A melhor opcao e: N=" + str(melhor_indice+1) + ", com r=" + repr(res[melhor_indice]))
'''

###### Gera previsao pros proximos 100 dias #######

print("")
print("######## Gera previsao ##########")
print("")

#descomentar, se comentar o codigo de calcular o melhor N
melhor_indice = N-1

#cria mais 100 particoes
particoesF = np.zeros((num_pontos,), dtype=np.float128)
particoesF = np.arange(b, b+100, dp)
particoesF = np.append(particoesF, b+100)

#calcula pontos interpolados

num_pontos_futuros = len(particoesF)

zf = np.zeros((num_pontos_futuros,), dtype=np.float128)

#print("num_pontos_futuros", num_pontos_futuros)

for k in range (0, num_pontos_futuros):
    #print(k)
    xk = particoesF[k]
    zf[k] = MinimosQuadrados().interpolaCoeficientes(coeficientes, melhor_indice+1, xk)
    
#print("coeficientes", coeficientes)
#print("particoesF", particoesF)
#print("zf", zf)

gerarGraficoPrevisao(x, y, particoes, z, particoesF, zf, melhor_indice+1)
#gerarGraficoPrevisao(x, y, particoes, z, [501, 502, 513, 544, 565], [26, 20, 15, 12, 14], N)

###### Dados futuros em relacao ao arquivo #######
#o arquivo para no dia 31/10/2018

#dados futuros (a partir de 01/11/2018)
xf = np.array([501, 502, 503, 504, 505, 506, 507], dtype=np.float128)
yf = np.array([27.62, 27.32, 27.32, 27.32, 27.32, 28.16, 27.19], dtype=np.float128)