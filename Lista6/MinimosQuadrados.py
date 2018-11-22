#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  3 21:06:23 2018

@author: Thiago Almeida
"""

import numpy as np
from metodos_numericos.LU import LU
from metodos_numericos.Gauss import Gauss
from Utils import Utils

class MinimosQuadrados():

    #constroi o polinomio para interpolacao
    def executar(self, a, b, n, f):
        
        #cria matriz
        for i in range (0,n):
            for j in range (0,n):
                fi_a = x**i
                fi_b = x**j
                
                A[i][j] = self.integral(a, b, fi(x, i), fi(x, j))
                
        #cria vetor fonte
        for i in range (0,n):
            b = x**i
            B[i] = self.integral(f, b)
            
        print("A", A)
        print("B", B)
        
        #Utils().obtemInfoMatriz(A)
            
        #calcula coeficientes
        X = LU().executar(A, B)[0]
        #X = Gauss().executarComPivoteamento(A, B)[0]
        #X = Gauss().executar(A, B)[0]
        return X
    def fi(x):
        return x**i
        
    def integral(self, a, b):
        return 0
    
    def interpolaCoeficientes(self, c, n, xk):
        
        soma = 0        
        for i in range (0,n+1):
            soma += c[i] * (xk ** i)
            
        #if(np.isnan(soma)):
            #print("resultado da interpolacao do valor "+repr(xk)+" foi igual a NaN")
            #soma = 0
            
        return soma
        