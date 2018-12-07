from metodos_numericos.RetroSubstituicao import RetroSubstituicao
from Utils import Utils
import numpy as np

class Gauss():
    def executar(self, M, B, dataType):
        if(len(M)==0):
            return 0

        self.dataType = dataType
        
        Ma = Utils().copiaMatriz(M)
        Ba = list(B)

        ordem = len(M[0])
        
        #percorre os pivos da matriz zerando as linhas abaixo
        for k in range(0, ordem):
            t = k + 1
            
            #percorre os elementos abaixo do pivo para extrair o multiplicador
            for i in range(t, ordem):  
                mult = Ma[i][k] / Ma[k][k]
                
                if(mult == 0 or np.isnan(mult)):
                    print("************** PIVO NULO!!!!!")
                if(mult >= 1):
                    print("MULTIPLICADOR MAIOR OU IGUAL A 1", mult)
                
                #usa o multiplicador pra zerar o elemento da linha
                for j in range(t, ordem):  
                    Ma[i][j] = Ma[i][j] - mult * Ma[k][j]
                
                #usa o multiplicador pra mudar o elemento no vetor fonte
                Ba[i] = Ba[i] - mult * Ba[k]
                
        #agora que tem uma matriz diagonal superior, usa retroSubstituicao
        retrosub = RetroSubstituicao().executar(Ma, Ba)
            
        return retrosub[0]

    def executarComPivoteamento(self, M, B, dataType):    
        ordem = len(M[0])
        
        self.dataType = dataType

        Ma = np.array(M, dtype=dataType, copy=True)#Utils().copiaMatriz(M)
        Ba = np.array(B, dtype=dataType, copy=True)#list(B)
        
        #loop do pivoteamento
        for a in range(0, ordem):

            #armazena pivo em modulo
            maior_elemento = abs(Ma[a][a])
            linha_a = a
            linha_b = a

            #percorre as colunas a procura do maior valor
            for b in range(a, ordem):
                elemento = abs(Ma[b][a])
                if(elemento > maior_elemento):
                    maior_elemento = elemento
                    linha_b = b
            
            #se os indices das linhas a e b forem diferentes, faz a troca de linha
            if(linha_a != linha_b):
                Utils().trocarLinhas(Ma, Ba, linha_a, linha_b)

            #percorre os elementos abaixo do pivo para extrair o multiplicador
            for i in range(a+1, ordem):  
                mult = Ma[i][a] / Ma[a][a]
                
                if(mult == 0 or np.isnan(mult)):
                    print("************** PIVO NULO!!!!!")
                if(mult >= 1):
                    print("MULTIPLICADOR MAIOR OU IGUAL A 1", mult)
                
                #usa o multiplicador pra zerar o elemento da linha
                for j in range(a, ordem):  
                    Ma[i][j] = Ma[i][j] - mult * Ma[a][j]
                    #if(passos % 100000 == 0):
                        #print("Passos percorridos: " + str(passos))
                
                #usa o multiplicador pra mudar o elemento no vetor fonte
                Ba[i] = Ba[i] - mult * Ba[a]

        #agora que tem uma matriz diagonal superior, usa retroSubstituicao
        retrosub = RetroSubstituicao().executar(Ma, Ba, dataType)
        return retrosub[0]