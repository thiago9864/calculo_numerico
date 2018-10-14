# -*- coding: utf-8 -*-

#from __future__ import division 
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import math as m
#import numpy as np
#from numpy import linalg  

#imports locais
from Utils import Utils
#from metodos_numericos.RetroSubstituicao import RetroSubstituicao
from metodos_numericos.Gauss import Gauss
from metodos_numericos.GaussSeidel import GaussSeidel
from metodos_numericos.Cholesky import Cholesky
from metodos_numericos.Thomas import Thomas
from metodos_numericos.Jacobi import Jacobi
from metodos_numericos.LU import LU

#talvez usar
#https://stackoverflow.com/questions/20548628/how-to-do-parallel-programming-in-python
#https://sebastianraschka.com/Articles/2014_multiprocessing.html

#Definicao da funcao solucao exata
def F(x, c1, c2, E): 
    return c1 * m.exp(- x / m.sqrt (E)) +  c2 * m.exp(x / m.sqrt (E)) + 1.0
       
##### codigo da questao 2 da lista 1 que gera a matriz do problema de valor de contorno #####

def gerarMatriz(num_particoes, E):

    #Coeficientes
    a = 0
    b = 1.0
    Nel = num_particoes
    ordem = num_particoes - 1
    h = (b-a)/float(Nel)
    c2 = (m.exp((-1.0)/m.sqrt(E)) - 1.0) / (m.exp(1.0/m.sqrt(E)) - m.exp((-1.0)/m.sqrt(E))) 
    c1 = - 1 - c2
    
    h2 = h**2.0
    cdp = 2 * E + h2

    
    #calcula tempo e solucao exata
    tempo = []
    solucao_exata = []
    
    for i in range(num_particoes+1):
        dh = h * i
        tempo.append(dh)
        solucao_exata.append(F(dh, c1, c2, E)) 


    #termo fonte
    B = []

    #monta matriz
    M = Utils().inicializaMatriz(ordem)

    for i in range(ordem):
        for j in range(ordem):
            if(i==j):
                #preenche diagonal principal
                M[i][j] = cdp
                B.append(h2)
            elif(i==j+1):
                #preenche diagonal inferior
                M[i][j] = E * -1.0
            elif(i==j-1):
                #preenche diagonal superior
                M[i][j] = E * -1.0
            else:
                #preenche o resto com zeros
                M[i][j] = 0
    
    return [M, B, tempo, solucao_exata]

            
##### Gerador de grafico #####

def gerarGrafico(tempo, solucao_aproximada, solucao_exata, metodo):
    
    #print("len tempo: ", len(tempo))
    #print("len solucao_aproximada: ", len(solucao_aproximada))
    #print("len solucao_exata: ", len(solucao_exata))
    
    #insere intervalos de contorno
    solucao_aproximada.insert(0, 0)
    solucao_aproximada.append(0)
    
     #plota grafico da função
    plt.plot(
        tempo, solucao_exata, 'b--',
        tempo, solucao_aproximada, 'r-'    
        )
    plt.ylabel(u"Valor de u(h)") #esse 'u' antes da string é pra converter o texto pra unicode
    plt.xlabel(u"Valor de h, " + str(len(tempo)-1) + u" partições")
    
    #legendas do grafico
    se_line = mlines.Line2D([], [], color='blue', marker='', markersize=0, label=u'Sol Exata')
    ac_line = mlines.Line2D([], [], color='red', marker='', markersize=0, label=u'Met. ' + metodo)
    
    plt.legend(handles=[se_line, ac_line])
    
    plt.title(u"Metodo de "+metodo+u" x Solução exata", )
    
    #plt.axis([0, 50, 0, 100])
    plt.show()   
    


##### Execucao dos codigos #####

numero_de_particoes = 1000
erro_do_metodo = 0.01

prev_passos = int((1.0/3.0) * (numero_de_particoes**3))

res = gerarMatriz(numero_de_particoes, erro_do_metodo)
M = res[0]
B = res[1]

#imprimeMatriz(M, B)

if(Utils().checarCriterioDasLinhas(M)):
    print("*** A matriz satisfaz o criterio das linhas ***")
else:
    print("*** A matriz não satisfaz o criterio das linhas ***")

print("")
print("------")
print("")

print("Metodo de Gauss (direto)")
inicio = Utils().getTime()
resGauss = Gauss().executar(M, B)
fim  = Utils().getTime()
print("Tamanho da matriz: " + repr(numero_de_particoes) + "x" + repr(numero_de_particoes))
print("Passos ate a resolucao: " + repr(resGauss[1]))
Utils().imprimeDiferencaTempo(inicio, fim)
gerarGrafico(res[2], resGauss[0], res[3], "Gauss")

print("")
print("------")
print("")

print("Metodo de Gauss Pivoteado Parcialmente (direto)")
print("previsao de passos: " + repr(prev_passos))
inicio = Utils().getTime()
resGauss = Gauss().executarComPivoteamento(M, B)
fim  = Utils().getTime()
print("Tamanho da matriz: " + repr(numero_de_particoes) + "x" + repr(numero_de_particoes))
print("Passos ate a resolucao: " + repr(resGauss[1]))
Utils().imprimeDiferencaTempo(inicio, fim)
gerarGrafico(res[2], resGauss[0], res[3], "Gauss Pivoteado Parcialmente")

print("")
print("------")
print("")

print("Metodo de Thomas (direto)")
inicio = Utils().getTime()
resThomas = Thomas().executar(M, B)
fim  = Utils().getTime()
print("Tamanho da matriz: " + repr(numero_de_particoes) + "x" + repr(numero_de_particoes))
print("Passos ate a resolucao: " + repr(resThomas[1]))
Utils().imprimeDiferencaTempo(inicio, fim)
gerarGrafico(res[2], resThomas[0], res[3], "Thomas")

print("")
print("------")
print("")

print("Metodo LU (direto)")
inicio = Utils().getTime()
resCholesky = LU().executar(M, B)
fim  = Utils().getTime()
print("Tamanho da matriz: " + repr(numero_de_particoes) + "x" + repr(numero_de_particoes))
print("Passos ate a resolucao: " + repr(resCholesky[1]))
Utils().imprimeDiferencaTempo(inicio, fim)
gerarGrafico(res[2], resCholesky[0], res[3], "LU")

print("")
print("------")
print("")
'''
print("Metodo de Cholesky (direto)")
inicio = Utils().getTime()
resCholesky = Cholesky().executar2(M, B)
fim  = Utils().getTime()
print("Tamanho da matriz: " + repr(numero_de_particoes) + "x" + repr(numero_de_particoes))
print("Passos ate a resolucao: " + repr(resCholesky[1]))
Utils().imprimeDiferencaTempo(inicio, fim)
gerarGrafico(res[2], resCholesky[0], res[3], "Cholesky")
'''
print("")
print("------")
print("")

print("Metodo de Jacobi (iterativo)")
chute_inicial = [1.0] * (numero_de_particoes - 1)
precisao = 0.00001
print("Erro esperado de: " + repr(precisao))
inicio = Utils().getTime()
resJacobi = Jacobi().executar2(M, B, chute_inicial, precisao, 5000)
fim  = Utils().getTime()
print("Tamanho da matriz: " + repr(numero_de_particoes) + "x" + repr(numero_de_particoes))
print("Iteracoes ate a resolucao: " + repr(resJacobi[2]))
print("Passos ate a resolucao: " + repr(resJacobi[1]))
Utils().imprimeDiferencaTempo(inicio, fim)
gerarGrafico(res[2], resJacobi[0], res[3], "Jacobi")

print("")
print("------")
print("")

print("Metodo de Gauss Seidel (iterativo)")
chute_inicial = [0] * (numero_de_particoes - 1)
precisao = 0.00001
inicio = Utils().getTime()
resGS = GaussSeidel().executar(M, B, chute_inicial, precisao, 1000)
fim  = Utils().getTime()
print("Tamanho da matriz: " + repr(numero_de_particoes) + "x" + repr(numero_de_particoes))
print("Passos ate a resolucao: " + repr(resGS[1]))
Utils().imprimeDiferencaTempo(inicio, fim)
gerarGrafico(res[2], resGS[0], res[3], "Gauss Seidel")
