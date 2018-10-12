from RetroSubstituicao import RetroSubstituicao
from Utils import Utils

class LU():
    def executar(self, M, B):
        if(len(M)==0):
            return 0

        ordem = len(M[0])
        passos = 0

        Ma = Utils().copiaMatriz(M)
        Ba = list(B)

        ##
        #aqui entra o codigo
        ##

        return [[0] * ordem, 0] 