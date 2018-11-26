#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  3 21:06:23 2018

@author: Thiago Almeida
"""

import numpy as np
from metodos_numericos.LU import LU
from metodos_numericos.Gauss import Gauss
#from Utils import Utils
#from TabelaGauss import TabelaGauss
from TabelaGaussLegendre import TabelaGaussLegendre

class MinimosQuadrados():

    #constroi o polinomio para interpolacao
    def executar(self, a, b, n, f):
        
        tam = n+1
        A = np.zeros((tam,tam), dtype=np.float128)
        B = np.zeros((tam,), dtype=np.float128)
        
        
        #cria matriz
        for i in range (0,tam):
            for j in range (0,tam):      
                A[i][j] = self.integralMatriz(a, b, tam, i, j)
                
        #cria vetor fonte
        for i in range (0,tam):
            B[i] = self.integralVetor(a, b, tam, f, i)
            
        #print("A", A)
        #print("B", B)
        
        #Utils().obtemInfoMatriz(A)
            
        #calcula coeficientes
        X = LU().executar(A, B)[0]
        #X = Gauss().executarComPivoteamento(A, B)[0]
        #X = Gauss().executar(A, B)[0]
        
        return X
    
    def funcaoBase(self, x, n):
        return x**n
    
    def produtoFi(self, x, i, j):
        return self.funcaoBase(x, i) * self.funcaoBase(x, j)
        
    def produtoVetor(self, x, f, i):
        return f(x) * self.funcaoBase(x, i)
    
    
    
    ######### Metodo de Integracao ##########
    
    
    
    
    def x(self, a, b, t):
        return (((b-a)*t) / 2.0) + ((b+a)/2.0)
        
    def dx(self, a, b):
        return (b-a)/2.0
        
    def integralMatriz(self, a, b, n, i, j):
        #recupera pontos da tabela
        tw = TabelaGaussLegendre().getValores(n)
        
        #calcula
        soma = 0
        for k in range (0, n):   
            wi = np.float128(tw[1][k])
            ti = np.float128(tw[0][k])

            soma += wi * self.produtoFi(self.x(a, b, ti), i, j) * self.dx(a, b)

        return soma
    
    def integralVetor(self, a, b, n, f, i):
        #recupera pontos da tabela
        tw = TabelaGaussLegendre().getValores(n)
        
        #calcula
        soma = 0
        for k in range (0, n):   
            wi = np.float128(tw[1][k])
            ti = np.float128(tw[0][k])

            soma += wi * self.produtoVetor(self.x(a, b, ti), f, i) * self.dx(a, b)

        return soma
    
    
    
    
    ######### Interpolacao ##########
    
    
    
    
    def interpolaCoeficientes(self, c, n, xk):
        
        tam = n+1
        
        soma = 0        
        for i in range (0,tam):
            soma += c[i] * (xk ** i)
            
        #if(np.isnan(soma)):
            #print("resultado da interpolacao do valor "+repr(xk)+" foi igual a NaN")
            #soma = 0
            
        return soma
        