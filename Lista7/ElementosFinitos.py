# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 11:22:02 2018

@author: thiagoalmeida
"""

import numpy as np
from TabelaGaussLegendre import TabelaGaussLegendre

class ElementosFinitos():
    
    '''
    Funcao base
    i = indice da funcao fi desejada
    n = numero de pontos    
    x = vetor com os pontos de x
    xk = ponto a avaliar
    '''
    def fi(self, i, n, x, xk):
        
        L = 0.0
        dL = 0.0
            
        #armazena os indices do vetor x pra usar na derivada
        ind = []
        #print("fi", i)
        #calcula o valor da funcao de base
        c = 1.0
        d = 1.0
        for j in range(n):
            if(i != j):
                ind.append(j)
                c = c * (xk - x[j])
                tmp = (x[i] - x[j])
                #print("tmp", tmp)
                d = d * tmp
        L = c / d
        
        #print("d", d)
        
        #calcula o valor da derivada da funcao de base
        soma = 0.0
        for j in range(n-1):
            k = ind[j]
            
            if(n > 2):
                soma = soma + (xk - x[k]) / d
            else:
                soma = soma + 1.0 / d
        
        dL = soma   
        
        #print("dL", dL)  

        return (L, dL)
    
   
    '''
    calcula o valor de K da matriz
    i = linha
    j = coluna
    h = tamanho do elemento
    pontos = numero de pontos do polinomio
    '''
    def K(self, i, j, h, pontos):
        
        print("i: " + str(i) + ", j: " + str(j))
        
        #recupera pontos da tabela
        tw = TabelaGaussLegendre().getValores(pontos)
        
        #gera pontos
        xp = np.zeros((pontos,), dtype=np.float128)
        
        q = 2.0 / (pontos-1)
        d = -1.0
        
        for k in range(pontos):
            xp[k] = d
            d += q
            
        #print("xp", xp)
        
        #calcula a integral por Gauss
        soma = 0
        for k in range (0, pontos):   
            ti = np.float128(tw[0][k])
            wi = np.float128(tw[1][k])
            
            #xk = self.x(a, b, ti)
            xk = ti
            
            (f1, df1) = self.fi(i, pontos, xp, xk)
            (f2, df2) = self.fi(j, pontos, xp, xk)
            

            f = df1 * (2.0/h) * df2 * (2.0/h)
            
            soma += wi * f * (h/2.0)
            
            print("ti", ti)
            print("wi", wi)
            print("df1", df1)
            print("df2", df2)
            print("h", h)
            
            print("-------")
            
        return soma
        
    
    '''
    calcula o valor de F do vetor fonte
    i = linha
    h = tamanho do elemento
    pontos = numero de pontos do polinomio
    '''
    def F(self, i, h, pontos):
        
        print("i: " + str(i))
        
        #recupera pontos da tabela
        tw = TabelaGaussLegendre().getValores(pontos)
        
        #gera pontos
        xp = np.zeros((pontos,), dtype=np.float128)
        
        q = 2.0 / (pontos-1)
        d = -1.0
        
        for k in range(pontos):
            xp[k] = d
            d += q
            
        #print("xp", xp)
        
        #calcula a integral por Gauss
        soma = 0
        for k in range (0, pontos):   
            ti = np.float128(tw[0][k])
            wi = np.float128(tw[1][k])
            
            #xk = self.x(a, b, ti)
            xk = ti
            
            (f1, df1) = self.fi(i, pontos, xp, xk)
            
            f = f1
            
            soma += wi * f * (h/2.0)
            
            print("ti", ti)
            print("wi", wi)
            print("f1", f1)
            print("h", h)
            
            print("-------")
            
        return soma
        
        
                    
            
            