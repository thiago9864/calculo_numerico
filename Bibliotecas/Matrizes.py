# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 14:17:33 2018

@author: thiagoalmeida
"""
import numpy as np
import math as m
import cmath as c

from Utils import Utils
from Datasource import Datasource


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
    Ma = Utils().copiaMatriz(M)
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



def gaussPivoteamento(M, B):    
    
    #cria uma copia
    Ma = Utils().copiaMatriz(M)
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
            Utils().trocarLinhas(Ma, Ba, linha_a, linha_b)
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
        erro = Utils().calculaErro(xp, x) 

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
        erro = Utils().calculaErro(X,Xa)
        
        if(erro < E):
            print("Terminou Gauss Seidel com erro de: ", erro)
            return [X, passos]
            
        #recebe vetor anterior
        Xa = list(X)#copia lista
    
    print("Gauss Seidel nao convergiu ou precisa de mais iteracoes para convergir")
    return [X, passos]
          
          
'''
def cholesky(A, B):
    
    #cria uma copia
    Ac = Utils().copiaMatriz(A)
    Bc = list(B)
    
    n = len(Ac[0])
    G = Utils().inicializaMatriz(n)
    
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
    
    #imprimeMatriz(G, [0] * n)
    
    Gt = Utils().copiaMatriz(G)
    G = Utils().transposicao(G)
    
    #imprimeMatriz(Gt, [0] * n)
    #imprimeMatriz(G, [0] * n)
    
    retrosub1 = retroSubstituicao(G, Bc)
    
    passos += retrosub1[1]
            
    #agora que tem uma matriz diagonal superior, usa retroSubstituicao
    retrosub2 = retroSubstituicao(Gt, retrosub1[0])
    return [retrosub2[0], retrosub2[1] + passos]
'''            
            
def cholesky(M, B):
    #cria uma copia
    A = Utils().copiaMatriz(M)
    B = np.array(B, np.float64)
    
    n = len(M)
    passos = 0
    
    L = Utils().inicializaMatriz(n)
    
    print("Matriz A")
    Utils().imprimeMatriz(A, B)
    
    ### Codigo copiado da Internet ###
    ### Fonte: https://rosettacode.org/wiki/Cholesky_decomposition#Python3.X_version_using_extra_Python_idioms
    L = [[0.0] * len(A) for _ in range(len(A))]
    for i in range(len(A)):
        for j in range(i+1):
            passos += 1
            s = sum(L[i][k] * L[j][k] for k in range(j))
            L[i][j] = m.sqrt(A[i][i] - s) if (i == j) else \
                      (1.0 / L[j][j] * (A[i][j] - s))

    Utils().imprimeMatriz(L, B)
    
    Lt = Utils().transposicao(L)
    
    Utils().imprimeMatriz(Lt, B)
          
    retrosub1 = retroSubstituicao(L, B)
    passos += retrosub1[1]
    
    #vetor resposta Y
    Y = retrosub1[0]

    retrosub2 = retroSubstituicao(Lt, Y)
    passos += retrosub1[1]
    
    #vetor resposta X
    X = retrosub2[0]
    
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





M = Datasource.m44C
B = Datasource.b44C

numero_de_particoes = len(M[0])


'''
inicio = Utils().getTime()
res = gauss(M, B)
fim  = Utils().getTime()
print("Resultado do sistema por Gauss: ", res[0])
print("Passos ate a resolucao: ", res[1])
Utils().imprimeDiferencaTempo(inicio, fim)
'''

print("--------------------")
inicio = Utils().getTime()
res = gaussPivoteamento(M, B)
fim  = Utils().getTime()
print("Resultado do sistema por Gauss (pivoteado): ", res[0])
print("Passos ate a resolucao: ", res[1])
Utils().imprimeDiferencaTempo(inicio, fim)


print("--------------------")
inicio = Utils().getTime()
res = cholesky(M, B)
fim  = Utils().getTime()
print("Resultado do sistema por Cholesky: ", res[0])
print("Passos ate a resolucao: ", res[1])
Utils().imprimeDiferencaTempo(inicio, fim)

'''
print("--------------------")
print("Metodo de Jacobi (iterativo)")
chute_inicial = [1.0] * numero_de_particoes
precisao = 0.01
inicio = Utils().getTime()
resJacobi = jacobi(M, B, chute_inicial, precisao, 100)
fim  = Utils().getTime()
print("Tamanho da matriz: " + repr(numero_de_particoes) + "x" + repr(numero_de_particoes))
print("Passos ate a resolucao: " + repr(resJacobi[1]))
Utils().imprimeDiferencaTempo(inicio, fim)
print("Resultado do sistema: ", resJacobi[0]);

print("--------------------")
print("Metodo de Gauss Seidel (iterativo)")
chute_inicial = [0] * numero_de_particoes
precisao = 0.01
inicio = Utils().getTime()
resGS = gaussSeidel(M, B, chute_inicial, precisao, 100)
fim  = Utils().getTime()
print("Tamanho da matriz: " + repr(numero_de_particoes) + "x" + repr(numero_de_particoes))
print("Passos ate a resolucao: " + repr(resGS[1]))
Utils().imprimeDiferencaTempo(inicio, fim)
print("Resultado do sistema: ", resGS[0])

''' 