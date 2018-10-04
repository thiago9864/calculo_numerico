# -*- coding: utf-8 -*-

import math as m
import sys
from datetime import datetime

#metodos uteis
def imprimeMatriz(matriz):
    ordem = len(matriz[0])
    for i in range(ordem):
        for j in range(ordem):
            sys.stdout.write(repr(matriz[i][j]).ljust(3))
        sys.stdout.flush()
        print("")

def inicializaMatriz(ordem):
    return [[0 for x in range(ordem)] for y in range(ordem)] 

def getTime():
    return datetime.now()
    
def imprimeDiferencaTempo(inicio, fim):
    dif = fim - inicio
    sys.stdout.write("Tempo de execucao: ")
    sys.stdout.flush()
    print(dif)

#metodos para calculo das condicoes dos metodos

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

##### codigo da questao 2 da lista 1 que gera a matriz do problema de valor de contorno #####


def gerarMatriz(num_particoes, E):

    #Coeficientes
    a = 0
    b = 1.0
    Nel = num_particoes
    h = (b-a)/float(Nel)
    c2 = (m.exp((-1.0)/m.sqrt(E)) - 1.0) / (m.exp(1.0/m.sqrt(E)) - m.exp((-1.0)/m.sqrt(E))) 
    c1 = - 1 - c2
    
    h2 = h**2.0
    cdp = 2 * E + h2

    B = []

    #monta matriz
    M = inicializaMatriz(num_particoes)

    for i in range(num_particoes):
        for j in range(num_particoes):
            if(i==j):
                #preenche diagonal principal
                M[i][j] = h2
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
            
            #usa o multiplicador pra mudar o elemento no vetor fonte
            B[i] = B[i] - mult * B[k]
            
    #agora que tem uma matriz diagonal superior, usa retroSubstituicao
    retrosub = retroSubstituicao(M, B)
    return [retrosub[0], retrosub[1] + passos]


##### Execucao dos codigos #####

ordem_da_matriz = 100
erro_do_metodo = 0.01

res = gerarMatriz(ordem_da_matriz, erro_do_metodo)
M = res[0]
B = res[1]

#imprimeMatriz(M)

inicio = getTime()
resGauss = gauss(M, B)
fim  = getTime()
#print("Resultado do sistema por Gauss: ", resGauss[0])
print("Tamanho da matriz: " + repr(ordem_da_matriz) + "x" + repr(ordem_da_matriz))
print("Passos ate a resolucao: " + repr(resGauss[1]))
imprimeDiferencaTempo(inicio, fim)