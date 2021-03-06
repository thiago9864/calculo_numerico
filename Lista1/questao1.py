# -*- coding: utf-8 -*-

import matplotlib.lines as mlines
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

def metodo_crank_nicolson(un, dt ,a):
    #return (2.0 / (2.0 + dt * a)) * ((2.0 * theta0 - dt * a * (un_ant - 2.0 * thetaM)) / 2.0)
    return (un * (2.0 - a * dt) + 2 * a * dt * thetaM) / (2.0 + a * dt)

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
        
    #metodo cranck-nicolson
    arr_crank_nicolson.append(vd)
    vd = metodo_crank_nicolson(vd, deltaT, K)     
        
    #arr_crank_nicolson.append(vd)
    #if c > 0:
        #calcular o proximo valor
    #    b = vd
    #    vd = metodo_crank_nicolson(deltaT, K, ub)
    #    #valor inicial = valor anterior
    #    ub = b
    #else:
    #    arr_crank_nicolson.append(vd)
    #    vd = vb
        
#mantem os arrays com mesmo tamanho reduzindo esse em 1
arr_dif_central.pop(1)
#arr_crank_nicolson.pop()
    
#print (arr_euler_explicito)

#plota grafico
plt.plot(
    arr_tempo, arr_solucao_exata, 'b--',
    arr_tempo, arr_euler_explicito, 'r-',
    arr_tempo, arr_euler_implicito, 'g-',
    arr_tempo, arr_dif_central, 'k-',
    arr_tempo, arr_crank_nicolson, 'c-'
    )
plt.ylabel(u"Temperatura (ºC)") #esse 'u' antes da string é pra converter o texto pra unicode
plt.xlabel(u"Tempo (seg), " + str(nParticoes) + u" partições")

#legendas do grafico
se_line = mlines.Line2D([], [], color='blue', marker='', markersize=0, label=u'Solução Exata')
ee_line = mlines.Line2D([], [], color='red', marker='', markersize=0, label=u'Euler Explícito')
ei_line = mlines.Line2D([], [], color='green', marker='', markersize=0, label=u'Euler Implícito')
dc_line = mlines.Line2D([], [], color='black', marker='', markersize=0, label=u'Diferença Central')
cn_line = mlines.Line2D([], [], color='cyan', marker='', markersize=0, label='Crank Nicolson')

plt.legend(handles=[se_line, ee_line, ei_line, dc_line, cn_line])

plt.title("Metodos de Resolucao")

#plt.axis([0, 50, 0, 100])
plt.show()
