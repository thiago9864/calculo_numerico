# -*- coding: utf-8 -*-

from metodos_numericos.RetroSubstituicao import RetroSubstituicao
from Utils import Utils
import math as m
import numpy as np

class Cholesky():
    #esse metodo esta com problema, pesquisar solucao
    #o que esta feito, esta baseado nesse link: https://math.stackexchange.com/questions/2422012/solving-a-linear-system-with-cholesky-factorization
    def executar(self, A, B):
        
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
        retrosub1 = RetroSubstituicao().executar(G, Bc)
        passos += retrosub1[1]
                
        #o resultado de cima e o vetor b do abaixo
        #Gtx=y
        retrosub2 = RetroSubstituicao().executar(Gt, retrosub1[0])
        return [retrosub2[0], retrosub2[1] + passos]
    

    
    def executar2(self, M, B):
        #cria uma copia
        A = Utils().copiaMatriz(M)
        B = np.array(B, np.float64)
        
        passos = 0
        
        ### Codigo copiado da Internet ###
        ### Fonte: https://rosettacode.org/wiki/Cholesky_decomposition#Python3.X_version_using_extra_Python_idioms
        ### Esse codigo retorna a matriz L
        L = [[0.0] * len(A) for _ in range(len(A))]
        for i in range(len(A)):
            for j in range(i+1):
                passos += 1
                s = sum(L[i][k] * L[j][k] for k in range(j))
                L[i][j] = m.sqrt(A[i][i] - s) if (i == j) else \
                          (1.0 / L[j][j] * (A[i][j] - s))
    
        #Resolucao do sistema baseado no link
        #https://math.stackexchange.com/questions/2422012/solving-a-linear-system-with-cholesky-factorization
        
        #transposta da matriz L
        Lt = Utils().transposicao(L)
        
        #resolve o primeiro passo do metodo
        retrosub1 = RetroSubstituicao().executar(L, B)
        passos += retrosub1[1]
        
        #vetor resposta Y
        Y = retrosub1[0]
    
        #resolve o segundo passo do metodo
        retrosub2 = RetroSubstituicao().executar(Lt, Y)
        passos += retrosub2[1]
        
        #vetor resposta X
        X = retrosub2[0]
        
        return [X, passos]

        ''' 
        Renan: Eu encontrei esse código aqui num arquivo da internet e fui comparando com o seu, tá igual, mas por via das dúvidas vou deixar implementado ele aqui em baixo
        #https://rosettacode.org/wiki/Cholesky_decomposition
        def cholesky(A):
        L = [[0.0] * len(M) for _ in range(len(M))]
        for i in range(len(M)):
            for j in range(i+1):
                s = soma(L[i][k] * L[j][k] for k in range(j))
                if (i == j):
                    L[i][j] = sqrt(M[i][i] - s) 
                else:1
                    L[i][j] = (1.0 / L[j][j] * (M[i][j] - s))
        '''

