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

#Matriz para cholesky (apostila)
m33C = [
    [ 4.0, -2.0,  2.0],
    [-2.0, 10.0, -7.0],
    [ 2.0, -7.0, 30.0]
]
b33C = [8.0, 11.0, -31.0]



#metodos para contagem de tempo
def getTime():
    return datetime.now()
    
def imprimeDiferencaTempo(inicio, fim):
    dif = fim - inicio
    sys.stdout.write("Tempo de execucao: ")
    sys.stdout.flush()
    print(dif)

def imprimeMatriz(matriz, vetor_solucao):
    ordem = len(matriz[0])
    just_space = 6
    for i in range(ordem):
        for j in range(ordem):
            sys.stdout.write("{0:.2f}".format(matriz[i][j]).ljust(just_space))
        sys.stdout.write("| " + repr(vetor_solucao[i]).ljust(just_space))
        sys.stdout.flush()
        print("")
    print("--------------------------")
    
def inicializaMatriz(ordem):
    return [[0 for x in range(ordem)] for y in range(ordem)] 
    
def copiaMatriz(matriz):
    Mr = list(matriz)
    ordem = len(matriz[0])
    for i in range(ordem):
        Mr[i] = list(matriz[i])
    return Mr
        
#Metodos necessarios para os algoritmos de resolucao de sistemas
def normaInfinito(x):
    size = len(x)
    maximo = abs(x[0])   
    
    for i in range(size):
        temp = abs(x[i])        
        if(temp > maximo):
                maximo = temp
            
    return maximo
    
def distanciaInfinito(x1, x2):
    if(len(x1) != len(x2)):
        print("O tamanho dos vetores x1 e x2 precisa ser o mesmo")
        return 0
        
    size = len(x1)
    dist = abs(x1[0] - x2[0])  
    
    
    for i in range(size):
        temp = abs(x1[0] - x2[0])
        if(dist > temp):
            temp = dist
            
    return dist
        
        
def calculaErro(x_prox, x_atual):
    return distanciaInfinito(x_prox, x_atual) / normaInfinito(x_prox)
    


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
    
    #cria uma copia
    Ma = copiaMatriz(M)
    Ba = list(B)
    
    if(len(Ma[0])==0):
        return 0
    
    ordem = len(Ma[0])
    passos = 0
    
    #percorre os pivos da matriz zerando as linhas abaixo
    for k in range(0, ordem):
        t = k + 1
        
        #percorre os elementos abaixo do pivo para extrair o multiplicador
        for i in range(t, ordem):  
            mult = Ma[i][k] / Ma[k][k]
            
            #usa o multiplicador pra zerar o elemento da linha
            for j in range(t, ordem):  
                Ma[i][j] = Ma[i][j] - mult * Ma[k][j]
                passos+=1
            
            #usa o multiplicador pra mudar o elemento no vetor fonte
            Ba[i] = Ba[i] - mult * Ba[k]
            
    #agora que tem uma matriz diagonal superior, usa retroSubstituicao
    retrosub = retroSubstituicao(Ma, Ba)
    return [retrosub[0], retrosub[1] + passos]

def trocarLinhas(Ma, Ba, linha_a, linha_b):
    
    ordem = len(Ma[0])
    #print("troca as linhas "+str(linha_a)+" e "+str(linha_b))
    for j in range(0, ordem):
        #troca a linha da matriz
        elemento = Ma[linha_a][j]
        Ma[linha_a][j] = Ma[linha_b][j]
        Ma[linha_b][j] = elemento

    #troca a linha do vetor fonte
    elemento = Ba[linha_a]
    Ba[linha_a] = Ba[linha_b]
    Ba[linha_b] = elemento

def gaussPivoteamento(M, B):    
    
    #cria uma copia
    Ma = copiaMatriz(M)
    Ba = list(B)
    
    ordem = len(Ma[0])
    passos = 0

    #imprimeMatriz(M, B)

    #loop do pivoteamento
    for a in range(0, ordem):

        #armazena pivo em modulo
        maior_elemento = abs(Ma[a][a])
        linha_a = a
        linha_b = a

        #percorre as colunas a procura do maior valor
        for b in range(a, ordem):
            passos += 1
            elemento = abs(Ma[b][a])
            if(elemento > maior_elemento):
                maior_elemento = elemento
                linha_b = b
        
        #se os indices das linhas a e b forem diferentes, faz a troca de linha
        if(linha_a != linha_b):
            trocarLinhas(Ma, Ba, linha_a, linha_b)
            #imprimeMatriz(M, B)

        #percorre os elementos abaixo do pivo para extrair o multiplicador
        for i in range(a+1, ordem):  
            mult = Ma[i][a] / Ma[a][a]
            
            #usa o multiplicador pra zerar o elemento da linha
            for j in range(a, ordem):  
                passos += 1
                Ma[i][j] = Ma[i][j] - mult * Ma[a][j]
            
            #usa o multiplicador pra mudar o elemento no vetor fonte
            Ba[i] = Ba[i] - mult * Ba[a]
        
        #print("calcula multiplicadores")
        #imprimeMatriz(M, B)

    #print("matriz triangular")
    #imprimeMatriz(M, B)

    #agora que tem uma matriz diagonal superior, usa retroSubstituicao
    retrosub = retroSubstituicao(Ma, Ba)
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
        x = list(xp)#copia a lista
            
    print("Jacobi nao convergiu ou precisa de mais iteracoes para convergir")
    return [xp, passos]

def gaussSeidel(M, B, chute_inicial, E, max_iteracoes):
    
    ordem = len(M[0])
    X = list(chute_inicial)
    Xa = list(chute_inicial)#vetor pra calcular o erro
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
        Xa = list(X)#copia lista
    
    print("Gauss Seidel nao convergiu ou precisa de mais iteracoes para convergir")
    return [X, passos]
          
          

def cholesky(A, B):
    
    #cria uma copia
    Ac = copiaMatriz(A)
    Bc = list(B)
    
    n = len(Ac[0])
    G = inicializaMatriz(n)
    
    passos = 0
    
    for i in range(n):
        soma = 0
        #gera a diagonal principal da matriz G
        for j in range(i+1):
            passos += 1
            if(j!=i):
                #soma os quadrados do que nao for da diagonal principal
                soma += G[i][j] ** 2
            else:
                #tira raiz quadrada da diferenca do item da diagonal principal de A
                #com a soma do que esta na G
                G[i][j] = m.sqrt(Ac[i][i] - soma) 
                
        #gera elementos fora da diagonal principal
        for i in range(j+1):
            soma = 0
            for k in range(j-1):
                passos += 1
                soma += G[i][k] * G[j][k]
            
            G[i][j] = (Ac[i][j] - soma) / G[j][j]
    
    #o resultado e uma matriz diagonal superior
    #que eu acho que e a G transposta
    #imprimeMatriz(G, [0] * n)
            
    #agora que tem uma matriz diagonal superior, usa retroSubstituicao
    retrosub = retroSubstituicao(G, Bc)
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



numero_de_particoes = len(m33J[0])

M = m33C;
B = b33C

imprimeMatriz(M, B)

inicio = getTime()
res = gauss(M, B)
fim  = getTime()
print("Resultado do sistema por Gauss: ", res[0]);
print("Passos ate a resolucao: ", res[1])
imprimeDiferencaTempo(inicio, fim)


print("--------------------")
inicio = getTime()
res = gaussPivoteamento(M, B)
fim  = getTime()
print("Resultado do sistema por Gauss (pivoteado): ", res[0]);
print("Passos ate a resolucao: ", res[1])
imprimeDiferencaTempo(inicio, fim)


print("--------------------")
inicio = getTime()
res = cholesky(M, B)
fim  = getTime()
print("Resultado do sistema por Cholesky: ", res[0]);
print("Passos ate a resolucao: ", res[1])
imprimeDiferencaTempo(inicio, fim)


print("--------------------")
print("Metodo de Jacobi (iterativo)")
chute_inicial = [1.0] * numero_de_particoes
precisao = 0.01
inicio = getTime()
resJacobi = jacobi(M, B, chute_inicial, precisao, 100)
fim  = getTime()
print("Tamanho da matriz: " + repr(numero_de_particoes) + "x" + repr(numero_de_particoes))
print("Passos ate a resolucao: " + repr(resJacobi[1]))
imprimeDiferencaTempo(inicio, fim)
print("Resultado do sistema: ", resJacobi[0]);

print("--------------------")
print("Metodo de Gauss Seidel (iterativo)")
chute_inicial = [0] * numero_de_particoes
precisao = 0.01
inicio = getTime()
resGS = gaussSeidel(M, B, chute_inicial, precisao, 100)
fim  = getTime()
print("Tamanho da matriz: " + repr(numero_de_particoes) + "x" + repr(numero_de_particoes))
print("Passos ate a resolucao: " + repr(resGS[1]))
imprimeDiferencaTempo(inicio, fim)
print("Resultado do sistema: ", resGS[0])

