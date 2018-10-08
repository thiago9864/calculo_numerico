# -*- coding: utf-8 -*-
from __future__ import division 
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import math as m
import numpy as np
from numpy import linalg  
import sys
from datetime import datetime

#talvez usar
#https://stackoverflow.com/questions/20548628/how-to-do-parallel-programming-in-python
#https://sebastianraschka.com/Articles/2014_multiprocessing.html

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
    
#Definicao da funcao solucao exata
def F(x, c1, c2, E): 
    return c1 * m.exp(- x / m.sqrt (E)) +  c2 * m.exp(x / m.sqrt (E)) + 1.0
    
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
    M = inicializaMatriz(ordem)

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
'''
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
'''

def jacobi(M, B, chute_inicial, E, max_iteracoes):

    
    x0 = chute_inicial
    nLinhas = len(M)
    x = 0.0 * len(x0)
    passos = 0
    
    if(verificaDiagonalDominante(M) == False):
        print("A matriz nao e diagonal dominante, portanto nao ira convergir")
        return [[0] * ordem, 0]
    
    for k in range(max_iteracoes):
        
        #interação do jaboci
        for i in range(nLinhas):
            x[i] = B[i]
            for j in range(1,i-1):
                passos+= 1 
                x[i] -= M[i][j] * x0[j]
            x[i] /= M[i][i]
        
        #se atingir o criterio de parada, interrompe e retorna os resultados
        erro = calculaErro(x0,x)
        
        if(erro < E):
            print("Terminou Jacobi com erro de: ", erro)
            return [x]
         
        #prepara proxima iteracao com aproximacao da anterior
        for i in range(len(x)):
            x0[i] = x[i]
            
    print("Jacobi nao convergiu ou precisa de mais iteracoes para convergir")
    
    return [x,passos]
    
    
    
    
    
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

numero_de_particoes = 10
erro_do_metodo = 0.01

#previsao para O(n^3)
prev_passos = int((2.0/3.0) * (numero_de_particoes**3))

res = gerarMatriz(numero_de_particoes, erro_do_metodo)
M = res[0]
B = res[1]

#imprimeMatriz(M)

'''
print("Metodo de Gauss (direto)")
inicio = getTime()
resGauss = gauss(M, B)
fim  = getTime()
print("Tamanho da matriz: " + repr(numero_de_particoes) + "x" + repr(numero_de_particoes))
print("Passos ate a resolucao: " + repr(resGauss[1]))
imprimeDiferencaTempo(inicio, fim)
gerarGrafico(res[2], resGauss[0], res[3], "Gauss")
'''

'''
print("Metodo de Gauss Pivoteado Parcialmente (direto)")
print("previsao de passos: " + repr(prev_passos))
inicio = getTime()
resGauss = gaussPivoteamento(M, B)
fim  = getTime()
print("Tamanho da matriz: " + repr(numero_de_particoes) + "x" + repr(numero_de_particoes))
print("Passos ate a resolucao: " + repr(resGauss[1]))
imprimeDiferencaTempo(inicio, fim)
gerarGrafico(res[2], resGauss[0], res[3], "Gauss Pivoteado Parcialmente")
'''

'''
print("Metodo de Thomas (direto)")
inicio = getTime()
resThomas = thomas(M, B)
fim  = getTime()
print("Tamanho da matriz: " + repr(numero_de_particoes) + "x" + repr(numero_de_particoes))
print("Passos ate a resolucao: " + repr(resThomas[1]))
imprimeDiferencaTempo(inicio, fim)
gerarGrafico(res[2], resThomas[0], res[3], "Thomas")
'''


print("Metodo de Jacobi (iterativo)")
chute_inicial = [0] * (numero_de_particoes - 1)
precisao = 0.01
inicio = getTime()
resJacobi = jacobi(M, B, chute_inicial, precisao, 1000)
fim  = getTime()
print("Tamanho da matriz: " + repr(numero_de_particoes) + "x" + repr(numero_de_particoes))
print("Passos ate a resolucao: " + repr(resJacobi[1]))
imprimeDiferencaTempo(inicio, fim)
gerarGrafico(res[2], resJacobi[0], res[3], "Jacobi")

