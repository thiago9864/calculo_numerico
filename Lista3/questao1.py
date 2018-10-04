# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import math as m
import sys
from datetime import datetime

#metodos uteis
def imprimeMatriz(matriz):
    ordem = len(matriz[0])
    for i in range(ordem):
        for j in range(ordem):
            sys.stdout.write(repr(matriz[i][j]).ljust(6))
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
    
#Definicao da funcao solucao exata
def F(x, c1, c2, E): 
    return c1 * m.exp(- x / m.sqrt (E)) +  c2 * m.exp(x / m.sqrt (E)) + 1.0
    

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
            
            #usa o multiplicador pra mudar o elemento no vetor fonte
            B[i] = B[i] - mult * B[k]
            
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

##### Gerador de grafico #####

def gerarGrafico(tempo, solucao_aproximada, solucao_exata, metodo):
    
    print("len tempo: ", len(tempo))
    print("len solucao_aproximada: ", len(solucao_aproximada))
    print("len solucao_exata: ", len(solucao_exata))
    
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

numero_de_particoes = 100
erro_do_metodo = 0.01

res = gerarMatriz(numero_de_particoes, erro_do_metodo)
M = res[0]
B = res[1]

#imprimeMatriz(M)

'''
print("Metodo de Gauss")
inicio = getTime()
resGauss = gauss(M, B)
fim  = getTime()
#print("Resultado do sistema por Gauss: ", resGauss[0])
print("Tamanho da matriz: " + repr(ordem_da_matriz) + "x" + repr(ordem_da_matriz))
print("Passos ate a resolucao: " + repr(resGauss[1]))
imprimeDiferencaTempo(inicio, fim)
gerarGrafico(res[2], resGauss[0], res[3], "Gauss")
'''

print("Metodo de Thomas")
inicio = getTime()
resThomas = thomas(M, B)
fim  = getTime()
#print("Resultado do sistema por Gauss: ", resGauss[0])
print("Tamanho da matriz: " + repr(ordem_da_matriz) + "x" + repr(ordem_da_matriz))
print("Passos ate a resolucao: " + repr(resThomas[1]))
imprimeDiferencaTempo(inicio, fim)
gerarGrafico(res[2], resThomas[0], res[3], "Thomas")
