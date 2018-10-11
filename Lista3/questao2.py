# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 11:29:03 2018

@author: ice
"""

# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import math as m
import sys
from datetime import datetime

#metodos uteis
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

def inicializaMatriz(ordem):
    return [[0 for x in range(ordem)] for y in range(ordem)] 

def getTime():
    return datetime.now()
    
def imprimeDiferencaTempo(inicio, fim):
    dif = fim - inicio
    sys.stdout.write("Tempo de execucao: ")
    sys.stdout.flush()
    print(dif)
    
    
#### Gera matriz para essa questao #####
    
    
def gerarMatriz(ordem):
    
    M = inicializaMatriz(ordem)
    B = [0] * ordem
    
    #cria a matriz
    for i in range(ordem):
        for j in range(ordem):
            M[i][j] = 1.0 / (i + j + 1.0)
    
    #cria o vetor fonte
    for i in range(ordem):
        B[i] = 1.0 / (i + ordem + 1.0)
        
    return [M, B]
   

##### Metodos de resolucao de sistemas #####


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
        div =  M[i][i]
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

        #percorre os elementos abaixo do pivo para extrair o multiplicador
        for i in range(a+1, ordem):  
            mult = M[i][a] / M[a][a]
            
            #usa o multiplicador pra zerar o elemento da linha
            for j in range(a, ordem):  
                passos += 1
                M[i][j] = M[i][j] - mult * M[a][j]
                #if(passos % 100000 == 0):
                    #print("Passos percorridos: " + str(passos))
            
            #usa o multiplicador pra mudar o elemento no vetor fonte
            B[i] = B[i] - mult * B[a]

    #agora que tem uma matriz diagonal superior, usa retroSubstituicao
    retrosub = retroSubstituicao(M, B)
    return [retrosub[0], retrosub[1] + passos]


def thomas(M, d):
    ''' Resolve Ax = d onde A e uma matriz tridiagonal composta pelos vetores a, b, c
		a - subdiagonal
        b - diagonal principal
		c - superdiagonal.
	Retorna x
    '''
    passos = 0
    num_particoes = len(M[0])
    a = [] #subdiagonal
    b = [] #diagonal principal
    c = [] #superdiagonal
    
    #separa os vetores da matriz
    for i in range(num_particoes):
        for j in range(num_particoes):
            passos += 1
            if(i==j):
                #preenche diagonal principal
                b.append(M[i][j])
            elif(i==j+1):
                #preenche diagonal inferior
                a.append(M[i][j])
            elif(i==j-1):
                #preenche diagonal superior
                c.append(M[i][j])


    n = len(d) # len(d) == len(b)
    c_ = [ c[0] / b[0] ]
    d_ = [ d[0] / b[0] ]
    
    for i in range(1, n):
        passos += 1
        aux = b[i] - c_[i-1]*a[i-1]
        if i < n-1:
            c_.append( c[i] / aux )
        d_.append( (d[i] - d_[i-1]*a[i-1])/aux )
    
    # Substituicao de volta
    x = [d_[-1]]
    for i in range(n-2, -1, -1):
        passos += 1
        x = [ d_[i] - c_[i]*x[0] ] + x
    
    return [x, passos]
    
    
##### Metodos Iterativos #####
    
def verificaDiagonalDominante(M):
    ordem = len(M)
    #percorre a matriz
    for i in range(ordem):
        soma = 0
        diag = M[i][i]
        for j in range(ordem):
            if(i != j):
                soma += abs(M[i][j])
        if(soma >= diag):
            return False
    return True

def jacobi(M, B, chute_inicial, E, max_iteracoes):
    ordem = len(M)
    x = chute_inicial
    xp = [0] * len(x)
    passos = 0
    
    if(verificaDiagonalDominante(M) == False):
        print("A matriz nao e diagonal dominante, portanto nao ira convergir")
        return [[0] * ordem, 0]
    
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
        for i in range(len(xp)):
            x[i] = xp[i]
            
    print("Jacobi nao convergiu ou precisa de mais iteracoes para convergir")
    return [xp, passos]
    
    
    
#### Executa os metodos #####

numero_de_elementos = 3  
    
mat = gerarMatriz(numero_de_elementos)

M = mat[0]
B = mat[1]

imprimeMatriz(M, B)

print("")
print("------")
print("")

print("Metodo de Gauss (direto)")
inicio = getTime()
resGauss = gauss(M, B)
fim  = getTime()
print("Tamanho da matriz: " + repr(numero_de_particoes) + "x" + repr(numero_de_particoes))
print("Passos ate a resolucao: " + repr(resGauss[1]))
imprimeDiferencaTempo(inicio, fim)
#print("Resultado do sistema por Gauss: ", resGauss[0]);


print("")
print("------")
print("")

print("Metodo de Jacobi (iterativo)")
chute_inicial = [0] * (numero_de_particoes - 1)
precisao = 0.01
inicio = getTime()
resJacobi = jacobi(M, B, chute_inicial, precisao, 1000)
fim  = getTime()
print("Tamanho da matriz: " + repr(numero_de_particoes) + "x" + repr(numero_de_particoes))
print("Passos ate a resolucao: " + repr(resJacobi[1]))
imprimeDiferencaTempo(inicio, fim)
#print("Resultado do sistema por Jacobi: ", resJacobi[0]);
