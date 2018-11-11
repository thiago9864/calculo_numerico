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
    def executar(self, x, y, n):
        
        n_col = len(x)
        n_linhas = n+1
        tam2 = n+1
        
        #ja preenche a primeira linha com 1
        M1 = np.ones((n_linhas,n_col), dtype=np.float128)

        for k in range (1,n_linhas):
            for j in range (0,n_col):
                M1[k][j] = x[j]**k
            
        #print("M1", M1)
        
        A = np.zeros((tam2,tam2), dtype=np.float128)
        B = np.zeros((tam2,), dtype=np.float128)
        
        #cria matriz
        for i in range (0,tam2):
            for j in range (0,tam2):
                a = np.array(M1[i], dtype=np.float128, copy=True)
                b = np.array(M1[j], dtype=np.float128, copy=True)
                
                A[i][j] = np.dot(a, b)
                
        #cria vetor fonte
        for i in range (0,tam2):
            b = np.array(M1[i], dtype=np.float128, copy=True)
            B[i] = np.dot(y, b)
            
        print("A", A)
        print("B", B)
        
        #Utils().obtemInfoMatriz(A)
            
        #calcula coeficientes
        #X = LU().executar(A, B)[0]
        X = Gauss().executarComPivoteamento(A, B)[0]
        #X = Gauss().executar(A, B)[0]
        return X
        
    
    
    def interpolaCoeficientes(self, c, n, xk):
        
        soma = 0        
        for i in range (0,n+1):
            soma += c[i] * (xk ** i)
            
        #if(np.isnan(soma)):
            #print("resultado da interpolacao do valor "+repr(xk)+" foi igual a NaN")
            #soma = 0
            
        return soma
        
    def calculaResiduo(self, y, x, n, c):
        
        tam = len(y)
        
        soma1 = 0
        soma2 = 0
        soma3 = 0
        for k in range (0, tam):
            soma1 += (y[k] - self.interpolaCoeficientes(c, n, x[k]))**2
            soma2 += y[k]**2
            soma3 += y[k]
        
        r2 = 1 - (soma1 / (soma2 - (1/tam)*soma3**2))
        
        if(r2 < 0):
            print("N="+repr(n)+", r^2 = " + repr(r2) + " e nao tem raiz real")
            return 0
        else:
            return np.sqrt(r2)            
