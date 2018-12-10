# -*- coding: utf-8 -*-

"""
@author: Renan, Thiago
"""

#Calculo do Erro

import numpy as np
#from MetodosIntegracao import NewtonCotes
#from TabelaGaussLegendre import TabelaGaussLegendre
#from MinimosQuadrados import MinimosQuadrados

class Erro():

    
    def normaMaximo(self, x):
        size = len(x)
        maximo = np.abs(x[0])   
        
        for i in range(size):
            temp = np.abs(x[i])        
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
        #print(type(x_prox[0]))
        #print(type(x_atual[0]))
        return self.distanciaMaximo(x_prox, x_atual) / self.normaMaximo(x_prox)

    
    '''
    Funcao de calculo de erro do slide
    nint: representa o numero de pontos de integracao de Gauss
    nen: numero de nos do elemento
    nel: numero de elementos
    exata: f(x) da solucao exata
    '''
    
    '''
    def w(self, l):
        return self.dataType(self.tw[1][l])
    
    def shg(self,i, j, k):
        return 0
    
    def erroL2_slide(self, nel, nen, nint, h, exata):
        erul2 = 0
        self.tw = TabelaGaussLegendre().getValores(nint)
        
        for n in range(1, nel):
            eru = 0
            for l in nint(1, nint):
                uh = 0
                xx = 0
                for i in range(1, nen):
                    uh = uh + self.shg(1, i, l) * u(i)
                    xx = xx + self.shg(1, i, l) * xl(i)
                    
                eru = eru + ((exata(xx) - uh)**2) * self.w(l) * h/2.0
                
            erul2 = erul2 + eru
            
        return np.sqrt(erul2)
    '''   
    
    
    
    
    