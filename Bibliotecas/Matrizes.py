# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 14:17:33 2018

@author: thiagoalmeida
"""
import numpy as np
import math as m
from datetime import datetime

#det(m22) tem q dar -2
m22 = [
    [1, 2],
    [3, 4]
]

#det(m33) tem q dar -10.10205
m33 = [
    [1,       0,       0],
    [-1.5,    0,       4.5],
    [1.12266, 2.2449, -3.3679]
]

#det(m55) tem q dar 148
m55 = [
    [1, 2, 5, 3, 2],
    [1, 3, 7, 3, 4],
    [0, 5, 2, 2, 1],
    [1, 3, 0, 1, 2],
    [0, 6, 7, 4, 7]
]

#Sistema da questão 1 da lista 2 de Mecânica de 2018.3
#respostas correspondem a [Fax, Fay, Faz, T, Fbx, Fbz]
#resposta: [-1.6250, -3.9542, -2.02, 3.6833, 0.5417, 5.2722]
m66 = [
    [1.0, 0.0, 0.0,  0.29412,  1.0,  0.0],
    [0.0, 1.0, 0.0,  0.36765,  0.0,  0.0],
    [0.0, 0.0, 1.0, -0.88235,  0.0,  1.0],
    [0.0, 0.0, 0.0, -2.20590,  0.0,  4.5],
    [0.0, 0.0, 0.0, -0.88233,  6.0,  0.0],
    [0.0, 0.0, 0.0, -1.10295, -4.5,  0.0]
]
b66 = [0, 2.6, 0, 15.6, 0, -6.5]


#Matriz diagonal superior para RetroSubstituição
m33Rs = [
    [1, 1, 1],
    [0, 2, 1],
    [0, 0, 3]
]
b33Rs = [3, 3, 3]

#Matriz diagonal inferior para RetroSubstituição
m44Ri = [
    [ 2.0,  0.0,  0.0,  0.0],
    [ 3.0,  5.0,  0.0,  0.0],    
    [ 1.0, -6.0,  8.0,  0.0],
    [-1.0,  4.0, -3.0,  9.0]
]
b44Ri = [4.0, 1.0, 48.0, 0.0]

#metodos para contagem de tempo
def getTime():
    return datetime.now()
    
def imprimeDiferencaTempo(inicio, fim):
    dif = fim - inicio
    print("Tempo ate solucao: ", dif)

def det(M):
    if(len(M)==0):
        return 0
        
    ordem = len(M[0])
    
    #Calcula determinante 2x2 simples em o(1)
    if(ordem == 2):
        return (M[0][0] * M[1][1]) - (M[0][1] * M[1][0])
    
    #Usa Sarry pra calcular o 3x3 em o(1)
    elif(ordem == 3):
        d1 = (M[0][0] * M[1][1] * M[2][2])
        d2 = (M[0][1] * M[1][2] * M[2][0])
        d3 = (M[0][2] * M[1][0] * M[2][1])
        
        d4 = (M[0][2] * M[1][1] * M[2][0])
        d5 = (M[0][0] * M[1][2] * M[2][1])
        d6 = (M[0][1] * M[1][0] * M[2][2])
        
        return (d1 + d2 + d3) - (d4 + d5 + d6) 
    
    #Usa metodo de Laplace matrizes nxn com n > 3 em 0(3)
    else:
        vdet = 0
        p = 1 #usa a linha 1
        
        #percorre a linha da matriz
        for n in range(1, ordem+1):            
            #separa a matriz auxiliar
            A = []
            for i in range(1, ordem+1):
                line = []
                for j in range(1, ordem+1):
                    
                    if(i != p and j != n):
                        line.append(M[i-1][j-1])
                if(len(line) > 0):
                    A.append(line)
                    
            #calcula o cofator
            cft = (-1.0)**(n+1) * M[p-1][n-1] * det(A)
            
            #concatena com o valor anterior
            vdet = vdet + cft
                
        return vdet
        
        
def retroSubstituicao(M, B):
    if(len(M)==0):
        return 0
    
    ordem = len(M[0])
    temp = 0
    passos = 0
    
    #cria array de solucao
    sol = [0] * ordem
    
    #checa se e superior
    isSuperior = (M[ordem-1][0] == 0)
    
    #percorre as linhas da matriz (ta funcionando)
    for n in range(0, ordem):  
        
        #define o i se for superior ou inferior
        if(isSuperior):
            i = ordem - n - 1
        else:
            i = n
        
        #valor do vetor solucao da linha correspondente
        temp = B[i]
        
        #soma os valores exceto o valor do pivo
        for j in range(0, ordem):
            if(j != i):
                temp += sol[j] * M[i][j] * -1
                passos += 1
                
       
        #interrompe se der divisao por zero
        div =  M[i][i];
        if(div == 0):
            return
            
        #calcula a solucao da linha, usando a soma e o valor do pivo
        sol[i] = temp / M[i][i]
        
    return [sol, passos]
        

def gauss(M, B):
    if(len(M)==0):
        return 0
    
    ordem = len(M[0])
    passos = 0
    
    #percorre os pivos da matriz zerando as linhas abaixo
    for k in range(0, ordem):
        t = k + 1
        
        #percorre os elementos abaixo do pivo para extrair o multiplicador
        for i in range(t, ordem):  
            mult = M[i][k] / M[k][k]
            
            #usa o multiplicador pra zerar o elemento da linha
            for j in range(t, ordem):  
                M[i][j] = M[i][j] - mult * M[k][j]
            
            #usa o multiplicador pra mudar o elemento no vetor fonte
            B[i] = B[i] - mult * B[k]
            
    #agora que tem uma matriz diagonal superior, usa retroSubstituicao
    retrosub = retroSubstituicao(M, B)
    return [retrosub[0], retrosub[1] + passos]
        

#vdet = det(m22)
#print("Determinante 2x2: " + str(vdet));

#vdet = det(m33)
#print("Determinante 3x3: " + str(vdet));

#vdet = det(m55)
#print("Determinante 5x5: " + str(vdet));

#res = retroSubstituicao(m33Rs, b33Rs)
#print("Resultado do sistema (diagonal superior): ", res[0]);
#print("Passos ate a resolucao: ", res[1])

#res = retroSubstituicao(m44Ri, b44Ri)
#print("Resultado do sistema (diagonal inferior): ", res[0]);
#print("Passos ate a resolucao: ", res[1])

inicio = getTime()
res = gauss(m66, b66)
fim  = getTime()
print("Resultado do sistema por Gauss: ", res[0]);
print("Passos ate a resolucao: ", res[1])
imprimeDiferencaTempo(inicio, fim)