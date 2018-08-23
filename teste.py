# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from math import *

#arrays
arr_solucao_exata = []
arr_euler_explicito = []
arr_euler_implicito = []
arr_tempo = []

#intervalo de tempo

nParticoes = 4**3
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
    return (1.0 - a * dt) * ini

def metodo_euler_implicito(ini, dt, a):
    return (1.0 / (1.0 + a * dt)) * ini

#main loop
va = vb = valorInicial

for c in range(nParticoes):

    t = float(deltaT) * c

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
    
print (arr_euler_explicito)

#plota grafico
plt.plot(
    arr_tempo, arr_solucao_exata, 'b-',
    arr_tempo, arr_euler_explicito, 'r-',
    arr_tempo, arr_euler_implicito, 'g-'
    )
plt.ylabel('teste')
plt.axis([0, 50, 0, 100])
plt.show()