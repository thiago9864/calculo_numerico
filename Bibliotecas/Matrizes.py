# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 14:17:33 2018

@author: thiagoalmeida
"""
import numpy as np
import math as m
import sys
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

#Matriz para pivoteamento
m33Pv = [
    [ 2.0,  4.0, -2.0],
    [ 4.0,  9.0, -3.0],
    [-2.0, -3.0,  7.0]
]
b33Pv = [2.0, 8.0, 10.0]

#Matriz para jacobi (apostila)
m33J = [
    [ 4.00,  0.24, -0.08],
    [ 0.09,  3.00, -0.15],
    [ 0.04, -0.08,  4.00]
]
b33J = [8.0, 9.0, 20.0]

#metodos para contagem de tempo
def getTime():
    return datetime.now()
    
def imprimeDiferencaTempo(inicio, fim):
    dif = fim - inicio
    print("Tempo ate solucao: ", dif)

def imprimeMatriz(matriz, vetor_solucao):
    ordem = len(matriz[0])
    just_space = 6
    for i in range(ordem):
        for j in range(ordem):
            sys.stdout.write(repr(matriz[i][j]).ljust(just_space))
        sys.stdout.write("| " + repr(vetor_solucao[i]).ljust(just_space))
        sys.stdout.flush()
        print("")
    print("--------------------------")

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
    
    ordem = len(M)
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
        div =  M[i][i]
        if(div == 0):
            return
            
        #calcula a solucao da linha, usando a soma e o valor do pivo
        sol[i] = temp / M[i][i]
        
    return [sol, passos]
        

def gauss(M, B):
    if(len(M)==0):
        return 0
    
    ordem = len(M)
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
                passos+=1
            
            #usa o multiplicador pra mudar o elemento no vetor fonte
            B[i] = B[i] - mult * B[k]
            
    #agora que tem uma matriz diagonal superior, usa retroSubstituicao
    retrosub = retroSubstituicao(M, B)
    return [retrosub[0], retrosub[1] + passos]

def trocarLinhas(M, B, linha_a, linha_b):
    ordem = len(M[0])
    #print("troca as linhas "+str(linha_a)+" e "+str(linha_b))
    for j in range(0, ordem):
        #troca a linha da matriz
        elemento = M[linha_a][j]
        M[linha_a][j] = M[linha_b][j]
        M[linha_b][j] = elemento

    #troca a linha do vetor fonte
    elemento = B[linha_a]
    B[linha_a] = B[linha_b]
    B[linha_b] = elemento

def gaussPivoteamento(M, B):    
    ordem = len(M[0])
    passos = 0

    #imprimeMatriz(M, B)

    #loop do pivoteamento
    for a in range(0, ordem):

        #armazena pivo em modulo
        maior_elemento = abs(M[a][a])
        linha_a = a
        linha_b = a

        #percorre as colunas a procura do maior valor
        for b in range(a, ordem):
            passos += 1
            elemento = abs(M[b][a])
            if(elemento > maior_elemento):
                maior_elemento = elemento
                linha_b = b
        
        #se os indices das linhas a e b forem diferentes, faz a troca de linha
        if(linha_a != linha_b):
            trocarLinhas(M, B, linha_a, linha_b)
            #imprimeMatriz(M, B)

        #percorre os elementos abaixo do pivo para extrair o multiplicador
        for i in range(a+1, ordem):  
            mult = M[i][a] / M[a][a]
            
            #usa o multiplicador pra zerar o elemento da linha
            for j in range(a, ordem):  
                passos += 1
                M[i][j] = M[i][j] - mult * M[a][j]
            
            #usa o multiplicador pra mudar o elemento no vetor fonte
            B[i] = B[i] - mult * B[a]
        
        #print("calcula multiplicadores")
        #imprimeMatriz(M, B)

    #print("matriz triangular")
    #imprimeMatriz(M, B)

    #agora que tem uma matriz diagonal superior, usa retroSubstituicao
    retrosub = retroSubstituicao(M, B)
    return [retrosub[0], retrosub[1] + passos]

##### Metodos Iterativos #####

def jacobi(M, B, chute_inicial, E, max_iteracoes):
    ordem = len(M[0])
    x = chute_inicial
    xp = [0] * len(x)
    passos = 0

    print("ordem da matriz", ordem)
    
    for k in range(max_iteracoes):
        
        #percorre a matriz
        for i in range(ordem):
            #comeca a soma pelo termo do vetor fonte
            soma = B[i]
            div = 0
            for j in range(ordem):
                passos += 1
                #separa o divisor
                if(i==j):
                    div = M[i][j]
                else:
                    soma += M[i][j] * x[j] * -1.0
            #cria vetor de solucoes para proxima iteracao com resultados da linha
            xp[i] = soma / div
        
        #se atingir o criterio de parada, interrompe e retorna os resultados
        erro = calculaErro(xp, x) 

        if(erro < E):
            print("Terminou Jacobi com erro de: ", erro)
            return [xp, passos]
        
        #prepara proxima iteracao com aproximacao da anterior
        x = xp.copy()
            
    print("Jacobi nao convergiu ou precisa de mais iteracoes para convergir")
    return [xp, passos]

def gaussSeidel(M, B, chute_inicial, E, max_iteracoes):
    
    x0 = chute_inicial
    ordem = len(M[0])
    X = chute_inicial.copy()
    Xa = chute_inicial.copy()#vetor pra calcular o erro
    passos = 0
    
    print("ordem da matriz", ordem)
    
    for k in range(max_iteracoes):
        
        #percorre a matriz
        for i in range(ordem):
            #comeca a soma pelo termo do vetor fonte
            soma = B[i]
            div = 0
            for j in range(ordem):
                passos += 1
                #separa o divisor
                if(i==j):
                    div = M[i][j]
                else:
                    soma += M[i][j] * X[j] * -1.0
            #cria vetor de solucoes para proxima iteracao com resultados da linha
            X[i] = soma / div
        
        #se atingir o criterio de parada, interrompe e retorna os resultados
        erro = calculaErro(X,Xa)
        
        if(erro < E):
            print("Terminou Gauss Seidel com erro de: ", erro)
            return [X, passos]
            
        #recebe vetor anterior
        Xa = X.copy()
    
    print("Gauss Seidel nao convergiu ou precisa de mais iteracoes para convergir")
    return [X, passos]
            
            
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

imprimeMatriz(m33J, b33J)

inicio = getTime()
res = gauss(m33J, b33J)
fim  = getTime()
print("Resultado do sistema por Gauss: ", res[0]);
print("Passos ate a resolucao: ", res[1])
imprimeDiferencaTempo(inicio, fim)
'''
inicio = getTime()
res = gaussPivoteamento(m33Pv, b33Pv)
fim  = getTime()
print("Resultado do sistema por Gauss (pivoteado): ", res[0]);
print("Passos ate a resolucao: ", res[1])
imprimeDiferencaTempo(inicio, fim)
'''

print("--------------------")
print("Metodo de Jacobi (iterativo)")
chute_inicial = [1.0] * len(m33J[0])
precisao = 0.01
inicio = getTime()
resJacobi = jacobi(m33J, b33J, chute_inicial, precisao, 100)
fim  = getTime()
print("Tamanho da matriz: " + repr(numero_de_particoes) + "x" + repr(numero_de_particoes))
print("Passos ate a resolucao: " + repr(resJacobi[1]))
imprimeDiferencaTempo(inicio, fim)
print("Resultado do sistema: ", resJacobi[0]);

print("--------------------")
print("Metodo de Gauss Seidel (iterativo)")
chute_inicial = [0] * len(m33J[0])
precisao = 0.01
inicio = getTime()
resGS = gaussSeidel(m33J, b33J, chute_inicial, precisao, 100)
fim  = getTime()
print("Tamanho da matriz: " + repr(numero_de_particoes) + "x" + repr(numero_de_particoes))
print("Passos ate a resolucao: " + repr(resGS[1]))
imprimeDiferencaTempo(inicio, fim)
print("Resultado do sistema: ", resGS[0])