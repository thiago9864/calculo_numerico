from RetroSubstituicao import RetroSubstituicao
from Utils import Utils

class Jacobi():
    def executar(self, M, B, chute_inicial, E, max_iteracoes):
        ordem = len(M)
        x = chute_inicial
        xp = [0] * len(x)
        passos = 0
        
        #if(Utils().verificaDiagonalDominante(M) == False):
            #print("A matriz nao e diagonal dominante, portanto nao ira convergir")
            #return [[0] * ordem, 0]
        
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