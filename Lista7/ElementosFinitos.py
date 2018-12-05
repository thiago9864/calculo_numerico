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
        
        #calcula o valor da funcao de base
        c = 1.0
        d = 1.0
        for j in range(n):
            if(i != j):
                ind.append(j)
                c = c * (xk - x[j])
                d = d * (x[i] - x[j])
        L = c / d
        
        #calcula o valor da derivada da funcao de base
        soma = 0.0
        for j in range(n-1):
            k = ind[j]
            soma = soma + (xk - x[k]) / d 
        
        dL = soma          

        return (L, dL)
    
    '''
    funcao q muda o intervalo da integral
    '''
    def x(self, a, b, t):
        return (((b-a)*t) / 2.0) + ((b+a)/2.0)
    
    '''
    calcula elemento diferencial
    '''
    def dx(self, a, b):
        return(b-a)/2.0
        
    '''
    calcula o valor de K da matriz
    i = linha
    j = coluna
    h = tamanho do elemento
    n = numero de pontos
    x = vetor com os pontos de x dos elementos
    '''
    def K(self, i, j, h, n, x):
        a = 0.0
        b = h
        ordem = n
        pontos = 2
        
        #recupera pontos da tabela
        tw = TabelaGaussLegendre().getValores(ordem)
        
        
        #calcula a integral por Gauss
        soma = 0
        for i in range (0, ordem):   
            wi = np.float128(tw[0][i])
            ti = np.float128(tw[1][i])
            
            xk = self.x(a, b, ti)
            
            (f1, df1) = self.fi(i, pontos, x, xk)
            (f2, df2) = self.fi(j, pontos, x, xk)
            
            f = df1 * df2
            
            soma += wi * f * self.dx(a, b)

        return soma
        
        
        
                    
            
            