from Utils import Utils

class GaussSeidel():

    def executar(self, M, B, chute_inicial, E, max_iteracoes):
        
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