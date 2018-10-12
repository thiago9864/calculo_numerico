# -*- coding: utf-8 -*-
from __future__ import division 
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import math as m
import numpy as np
from numpy import linalg  
from metodos_numericos.RetroSubstituicao import RetroSubstituicao

#imports locais
from Utils import Utils

#talvez usar
#https://stackoverflow.com/questions/20548628/how-to-do-parallel-programming-in-python
#https://sebastianraschka.com/Articles/2014_multiprocessing.html

#Definicao da funcao solucao exata
def F(x, c1, c2, E): 
    return c1 * m.exp(- x / m.sqrt (E)) +  c2 * m.exp(x / m.sqrt (E)) + 1.0
       
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
    M = Utils().inicializaMatriz(ordem)

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
    retrosub = RetroSubstituicao().executar(M, B)
        
    return [retrosub[0], retrosub[1] + passos]


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
    retrosub = RetroSubstituicao().retroSubstituicao(M, B)
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
    
#esse metodo esta com problema, pesquisar solucao
#o que esta feito, esta baseado nesse link: https://math.stackexchange.com/questions/2422012/solving-a-linear-system-with-cholesky-factorization
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
    
    Gt = Utils().copiaMatriz(G)
    G = Utils().transposicao(G)#transposta da transposta resulta na matriz original
    
    #imprimeMatriz(Gt, [0] * n)
    #imprimeMatriz(G, [0] * n)
    
    # Gy=b 
    retrosub1 = retroSubstituicao(G, Bc)
    passos += retrosub1[1]
            
    #o resultado de cima e o vetor b do abaixo
    #Gtx=y
    retrosub2 = retroSubstituicao(Gt, retrosub1[0])
    return [retrosub2[0], retrosub2[1] + passos]
    
    
##### Metodos Iterativos #####
    

def jacobi(M, B, chute_inicial, E, max_iteracoes):
    ordem = len(M)
    x = chute_inicial
    xp = [0] * len(x)
    passos = 0
    
    if(Utils().verificaDiagonalDominante(M) == False):
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
        erro = Utils().calculaErro(xp, x) 

        if(erro < E):
            print("Terminou Jacobi com erro de: ", erro)
            return [xp, passos]
        
        #prepara proxima iteracao com aproximacao da anterior
        x = list(xp)
            
    print("Jacobi nao convergiu ou precisa de mais iteracoes para convergir")
    return [xp, passos]


#Thiago (funcionou, mas esta dando muitos passos)
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
#Renan
def gaussSeidel(M,B,chute_inicial,E):
    x0 = chute_inicial
    passos = 0
    N = len(M)
    
    #inicia o vetor
    for k in range(len(B)):
        X[k] = 0
    
    k = 1
    while 1 <= N:
        
        for i in range(N):
            passos += 1
            soma = 0
            
            for j in range(i-1):
                soma = soma + M[i][j] * x[j]
                
            for j in range(i+1,N):
                soma = soma + M[i][j] * X[j]
                
        x[i] = (b[i] - soma)/M[i][i]
        
        erro = calculaErro(X,x)
        
        if(erro < E):
            print("Terminou Gauss-Seidel com erro de:", erro)
            return [x]
        
        X[i] = x[i]
        
        return [x,passos]
'''
        
            
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

numero_de_particoes = 20
erro_do_metodo = 0.01

#previsao para O(n^3)
prev_passos = int((2.0/3.0) * (numero_de_particoes**3))

res = gerarMatriz(numero_de_particoes, erro_do_metodo)
M = res[0]
B = res[1]

#imprimeMatriz(M, B)


print("Metodo de Gauss (direto)")
inicio = Utils().getTime()
resGauss = gauss(M, B)
fim  = Utils().getTime()
print("Tamanho da matriz: " + repr(numero_de_particoes) + "x" + repr(numero_de_particoes))
print("Passos ate a resolucao: " + repr(resGauss[1]))
Utils().imprimeDiferencaTempo(inicio, fim)
gerarGrafico(res[2], resGauss[0], res[3], "Gauss")


'''
print("Metodo de Gauss Pivoteado Parcialmente (direto)")
print("previsao de passos: " + repr(prev_passos))
inicio = Utils().getTime()
resGauss = gaussPivoteamento(M, B)
fim  = Utils().getTime()
print("Tamanho da matriz: " + repr(numero_de_particoes) + "x" + repr(numero_de_particoes))
print("Passos ate a resolucao: " + repr(resGauss[1]))
Utils().imprimeDiferencaTempo(inicio, fim)
gerarGrafico(res[2], resGauss[0], res[3], "Gauss Pivoteado Parcialmente")
'''

'''
print("Metodo de Thomas (direto)")
inicio = Utils().getTime()
resThomas = thomas(M, B)
fim  = Utils().getTime()
print("Tamanho da matriz: " + repr(numero_de_particoes) + "x" + repr(numero_de_particoes))
print("Passos ate a resolucao: " + repr(resThomas[1]))
Utils().imprimeDiferencaTempo(inicio, fim)
gerarGrafico(res[2], resThomas[0], res[3], "Thomas")
'''

'''
print("Metodo de Cholesky (direto)")
inicio = Utils().getTime()
resCholesky = cholesky(M, B)
fim  = Utils().getTime()
print("Tamanho da matriz: " + repr(numero_de_particoes) + "x" + repr(numero_de_particoes))
print("Passos ate a resolucao: " + repr(resCholesky[1]))
Utils().imprimeDiferencaTempo(inicio, fim)
gerarGrafico(res[2], resCholesky[0], res[3], "Cholesky")
'''
'''
print("Metodo de Jacobi (iterativo)")
chute_inicial = [0.8] * (numero_de_particoes - 1)
precisao = 0.01
inicio = Utils().getTime()
resJacobi = jacobi(M, B, chute_inicial, precisao, 1000)
fim  = Utils().getTime()
print("Tamanho da matriz: " + repr(numero_de_particoes) + "x" + repr(numero_de_particoes))
print("Passos ate a resolucao: " + repr(resJacobi[1]))
Utils().imprimeDiferencaTempo(inicio, fim)
gerarGrafico(res[2], resJacobi[0], res[3], "Jacobi")
'''
'''
print("Metodo de Gauss Seidel (iterativo)")
chute_inicial = [0] * (numero_de_particoes - 1)
precisao = 0.01
inicio = Utils().getTime()
resGS = gaussSeidel(M, B, chute_inicial, precisao, 1000)
fim  = Utils().getTime()
print("Tamanho da matriz: " + repr(numero_de_particoes) + "x" + repr(numero_de_particoes))
print("Passos ate a resolucao: " + repr(resGS[1]))
Utils().imprimeDiferencaTempo(inicio, fim)
gerarGrafico(res[2], resGS[0], res[3], "Gauss Seidel")
'''