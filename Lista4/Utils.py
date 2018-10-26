# -*- coding: utf-8 -*-

import sys
from datetime import datetime

class Utils():
    def getTime(self):
        return datetime.now()
    
    def imprimeDiferencaTempo(self, inicio, fim):
        dif = fim - inicio
        sys.stdout.write("Tempo de execucao: ")
        sys.stdout.flush()
        print(dif)

    def imprimeMatriz(self, matriz):
            ordem = len(matriz[0])
            just_space = 8
            for i in range(ordem):
                for j in range(ordem):
                    sys.stdout.write("{0:.4f}".format(matriz[i][j]).ljust(just_space))
                sys.stdout.flush()
                print("")
            print("--------------------------")
                
    ######## Calculo do erro dos metodos numericos ########

    def normaMaximo(self, x):
        size = len(x)
        maximo = abs(x[0])   
        
        for i in range(size):
            temp = abs(x[i])        
            if(temp > maximo):
                maximo = temp
                
        return maximo
        
    def distanciaMaximo(self, x1, x2):
        if(len(x1) != len(x2)):
            print("O tamanho dos vetores x1 e x2 precisa ser o mesmo")
            return 0
            
        size = len(x1)
        dist = abs(x1[0] - x2[0])  
        
        for i in range(size):
            temp = abs(x1[i] - x2[i])
            if(temp > dist):
                dist = temp
                
        return dist
            
    def calculaErro(self, x_prox, x_atual):
        return self.distanciaMaximo(x_prox, x_atual) / self.normaMaximo(x_prox)

    def erroResidual(self, M, X, B):
        size = len(M[0])
        erroRet = []
        for i in range(size):
            valor = 0
            for j in range(size):
                valor += M[i][j] * X[j]

            erroRet.append(abs(valor - B[i]))

        return [erroRet, max(erroRet)]
