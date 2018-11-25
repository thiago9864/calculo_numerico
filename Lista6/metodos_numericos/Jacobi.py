#from metodos_numericos.RetroSubstituicao import RetroSubstituicao
from Utils import Utils
from numpy import zeros, diag, diagflat, dot
import numpy as np

class Jacobi():
    def executar(self, M, B, chute_inicial, E, max_iteracoes):
        ordem = len(M)
        x = chute_inicial
        xp = [0] * len(x)
        passos = 0
        iteracoes = 0
        
        #if(Utils().verificaDiagonalDominante(M) == False):
            #print("A matriz nao e diagonal dominante, portanto nao ira convergir")
            #return [[0] * ordem, 0]
        
        for k in range(max_iteracoes):
            iteracoes += 1
            
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
        return [xp, passos, iteracoes]


    def executar2(self, M, B, chute_inicial, E, max_iteracoes):
        ordem = len(M)
        X = np.array(chute_inicial, np.float64)
        Xp = np.array(chute_inicial, np.float64)
        B = np.array(B, np.float64)
        passos = 0
        iteracoes = 0
        
        for k in range(max_iteracoes):
            iteracoes += 1
            
            for i in range(ordem):
                soma = 0.0
                passos += 1
                #passa linha pra um array
                L = np.array(M[i], np.float64)
                #zera o lugar onde seria o valor da diagonal principal
                L[i] = 0
                #faz produto escalar entre os vetores
                soma = dot(L, X)
                '''
                for j in range(ordem):
                    passos += 1
                    if(i != j):
                        m = M[i][j]
                        x = X[j]
                        if(m != 0 and x != 0):
                            soma += (m * x)
                '''
                Xp[i] = (1.0 / M[i][i]) * (B[i] - soma)
             
            #se atingir o criterio de parada, interrompe e retorna os resultados
            erro = Utils().calculaErro(Xp, X) 

            if(erro < E):
                print("Terminou Jacobi com erro de: ", erro)
                return [list(Xp), passos, iteracoes]
            
            #prepara proxima iteracao com aproximacao da anterior
            X = np.array(Xp, np.float64)
            
            
            if(k % 100 == 0):
                print("Jacobi fez " + str(k) + " iteracoes...")
        
            
        print("Jacobi nao convergiu ou precisa de mais iteracoes para convergir")
        return [list(Xp), passos, iteracoes]
    
                
    def executar3(self, A,b,x,erro,N):
        """Solves the equation Ax=b via the Jacobi iterative method."""
        iteracoes = 0
        
        # Create an initial guess if needed                                                                                                                                                            
        if x is None:
            x = zeros(len(A[0]))
    
        # Create a vector of the diagonal elements of A                                                                                                                                                
        # and subtract them from A                                                                                                                                                                     
        D = diag(A)
        R = A - diagflat(D)
    
        # Iterate for N times                                                                                                                                                                          
        for i in range(N):
            iteracoes += 1
            x = (b - dot(R,x)) / D
        return [list(x), 0, iteracoes]

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