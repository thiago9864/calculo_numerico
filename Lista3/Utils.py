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


    ######## Metodos para manipulacao de matrizes ########


    def imprimeMatriz(self, matriz, vetor_solucao):
        ordem = len(matriz[0])
        just_space = 8
        for i in range(ordem):
            for j in range(ordem):
                sys.stdout.write("{0:.4f}".format(matriz[i][j]).ljust(just_space))
            sys.stdout.write("| " + repr(vetor_solucao[i]).ljust(just_space))
            sys.stdout.flush()
            print("")
        print("--------------------------")

    def inicializaMatriz(self, ordem):
        return [[0 for x in range(ordem)] for y in range(ordem)] 
        
    def copiaMatriz(self, matriz):
        Mr = list(matriz)
        ordem = len(matriz[0])
        for i in range(ordem):
            Mr[i] = list(matriz[i])
        return Mr
        
    def transposicao(self, matriz):
        ordem = len(matriz[0])
        Mt = self.inicializaMatriz(ordem)    
        for i in range(ordem):
            for j in range(ordem):
                Mt[i][j] = matriz[j][i]   
        return Mt

    def zeraMatriz(self, A):
        n = len(A[0])
        for i in range(n):
            for j in range(n):
                A[i][j] = 0

    def trocarLinhas(self, Ma, Ba, linha_a, linha_b):
    
        ordem = len(Ma[0])
        for j in range(0, ordem):
            #troca a linha da matriz
            elemento = Ma[linha_a][j]
            Ma[linha_a][j] = Ma[linha_b][j]
            Ma[linha_b][j] = elemento

        #troca a linha do vetor fonte
        elemento = Ba[linha_a]
        Ba[linha_a] = Ba[linha_b]
        Ba[linha_b] = elemento
                
    ######## Calculo do erro dos metodos numericos ########

    def normaInfinito(self, x):
        size = len(x)
        maximo = abs(x[0])   
        
        for i in range(size):
            temp = abs(x[i])        
            if(temp > maximo):
                    maximo = temp
                
        return maximo
        
    def distanciaInfinito(self, x1, x2):
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
            
    def calculaErro(self, x_prox, x_atual):
        return self.distanciaInfinito(x_prox, x_atual) / self.normaInfinito(x_prox)

    def erroResidual(self, M, X, B):
        size = len(M[0])
        erroRet = []
        for i in range(size):
            valor = 0
            for j in range(size):
                valor += M[i][j] * X[j]

            erroRet.append(abs(valor - B[i]))

        return [erroRet, max(erroRet)]

    ######## Metodos para verificacao de matrizes ########

    def det(self, M):
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
                cft = (-1.0)**(n+1) * M[p-1][n-1] * self.det(A)
                
                #concatena com o valor anterior
                vdet = vdet + cft
                    
            return vdet

    def verificaDiagonalDominante(self, M):
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

    def verificaPositivaDefinida(self, M):
        ordem = len(M)
        
        #define as k submatrizes
        for k in range(1, ordem):
            
            #inicializa a  submatriz
            A = self.inicializaMatriz(k)
            
            #monta a submatriz
            for i in range(k):
                for j in range(k):
                    A[i][j] = M[i][j]
                    
            #verifica se o det(A) <= 0, se for retorna false
            if(self.det(M) <= 0):
                return False
                
        #se passar tudo, e positiva definida
        return True


    def obtemInfoMatriz(self, M):
        det = self.det(M)
        print("----")
        if(det > 0):
            print("Determinante da matriz é maior que zero, e igual a: " + repr(det))
        else:
            print("Determinante da matriz é menor que zero, e igual a: " + repr(det))

        if(self.verificaDiagonalDominante(M)):
            print("A matriz é diagonal dominante")
        else:
            print("A matriz não é diagonal dominante")

        if(self.verificaPositivaDefinida(M)):
            print("A matriz é positiva definida")
        else:
            print("A matriz não é positiva definida")

    def checarCriterioDasLinhas(self, M):

        ordem = len(M)

        for i in range(ordem):
            valores = [] 
            div = M[i][i]

            #se algum elemento da diagonal principal for zero
            #a matriz nao satisfaz o criterio das linhas
            if(div == 0):
                return False

            for j in range(ordem):
                if(i != j):
                    valores.append(M[i][j] / div)
                
            #um elemento dividido pelo valor da diagonal principal deu maior ou igual que 1
            #a matriz nao satisfaz o criterio das linhas
            if(max(valores) >= 1):
                return False

        return True