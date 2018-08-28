# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from math import *

#arrays
arr_solucao_exata = []
arr_euler_explicito = []
arr_euler_implicito = []
arr_dif_central = []
arr_crank_nicolson = []
arr_tempo = []

#intervalo de tempo

nParticoes = 4**6
rangeT = range(0, 50)

#coeficientes
K = 0.035871952
theta0 = 99.0
thetaM = 27.0
valorInicial = theta0
deltaT = len(rangeT) / float(nParticoes)

print("Numero de Partições = " + str(nParticoes))
print("DeltaT = " + str(deltaT))
#verifica deltaT
if deltaT < (1 / K):
    print("DeltaT está ok")
else:
    print("DeltaT fora da especificacao")

#metodos para resolucao
def solucao_exata(t):
    return (theta0 - thetaM) * exp(K * t * -1) + thetaM

#metodos de aproximação
def metodo_euler_explicito(ini, dt, a):
    return (-1 * a * ini + a * thetaM) * dt + ini

def metodo_euler_implicito(ini, dt, a):
    return ((a * dt * thetaM) + ini) / (1 + a * dt)
    
def metodo_diferenca_central(un, un_ant, dt, a):
    return 2 * a * dt * (thetaM - un) + un_ant

#main loop
va = vb = vc = vd = valorInicial
ua = valorInicial

for c in range(nParticoes):

    t = deltaT * c

    #calcula solucao exata
    sol = solucao_exata(t)
    arr_solucao_exata.append(sol)
    arr_tempo.append(t)

    #metodo euler explicito
    arr_euler_explicito.append(va)
    va = metodo_euler_explicito(va, deltaT, K)

    #metodo euler implicito
    arr_euler_implicito.append(vb)
    vb = metodo_euler_implicito(vb, deltaT, K)
    
    #metodo diferenca central
    arr_dif_central.append(vc)
    if c > 0:
        #calcula proximo valor
        a = vc
        vc = metodo_diferenca_central(vc, ua, deltaT, K)
        #valor inicial = valor anterior
        ua = a
    else:
        arr_dif_central.append(vb)
        vc = vb
        
        
#mantem os arrays com mesmo tamanho reduzindo esse em 1
arr_dif_central.pop()
    
#print (arr_euler_explicito)

#plota grafico
plt.plot(
    arr_tempo, arr_solucao_exata, 'b--',
    arr_tempo, arr_euler_explicito, 'r-',
    arr_tempo, arr_euler_implicito, 'g-',
    arr_tempo, arr_dif_central, 'k-'
    )
plt.ylabel(u"Temperatura (ºC)") #esse 'u' antes da string é pra converter o texto pra unicode
plt.xlabel(u"Tempo (seg), " + str(nParticoes) + u" partições")
#plt.axis([0, 50, 0, 100])#isso aqui trava o gráfico em um retangulo com as proporções dadas em [l, t, r, b]
plt.show()