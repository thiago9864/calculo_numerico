# -*- coding: utf-8 -*-

#from metodos_numericos.RetroSubstituicao import RetroSubstituicao
from Utils import Utils
import numpy as np

class LU():
    
    def decomposicao(self,M):
        
        ordem = len(M[0])
        passos = 0
        
        #L = [[0.0 for i in range(ordem)]for j in range(ordem)]
        #U = [[0.0 for i in range(ordem)]for j in range(ordem)]
        
        L = Utils().iniciaMatrizComDataType(ordem, self.dataType)
        U = Utils().iniciaMatrizComDataType(ordem, self.dataType)
        
        #print(type(L[0][0]))
        #print(type(U[0][0]))
        
        for j in range(ordem):
            U[0][j] = M[0][j]
        for i in range(ordem):
            L[i][0] = M[i][0]/U[0][0]
        
        for i in range(ordem):
        #Calcula L
            for j in range(i+1):
                soma = 0.0
                for k in range(j):
                    passos+=1
                    soma += L[i][k] * U[k][j]
                L[i][j] = M[i][j] - soma
                    
        #Calcula U
            for j in range(i,ordem):
                soma = 0.0
                for k in range(i):
                   passos+=1
                   soma += L[i][k] * U[k][j]
                U[i][j] = (M[i][j] - soma)/L[i][i]
        return (L, U, passos)

    
    def executar(self, M, B, dataType):
        if(len(M)==0):
            return 0
        
        self.dataType = dataType
        
        ordem = len(M[0])
        
        y = np.zeros((ordem,), dtype=self.dataType)
        x = np.zeros((ordem,), dtype=self.dataType)
        (L, U, passos) = self.decomposicao(M)
        
        #print("L", L)
        #print("U", U)
        
        #Passo 1 para a resolução L * y = b
        
        y[0] = B[0]/L[0][0]
        for i in range(1, ordem):
            soma = 0.0
            for j in range(i):
                soma += L[i][j] * y[j]
            y[i] = (B[i] - soma)/L[i][i]
            
        #print("Y", y)
        
        #Passo 2 para a resolucao U * x = y
        x[ordem-1] = y[ordem-1]/U[ordem-1][ordem-1]
        for i in range(ordem-1, -1, -1):
            soma = y[i]
            for j in range(i+1, ordem):
                soma = soma - U[i][j] * x[j]
            x[i] = soma/U[i][i]
            
        #print("resultado lu:",x)
        
        return x   