#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 21:09:49 2018

@author: Thiago Almeida
"""
import numpy as np
#from Utils import Utils

class Newton():
    
    #cria uma variavel global
    coeficientes = 0
    
    def diferenca(self, xi, xf, fxi, fxf):
        return (fxf - fxi) / (xf - xi)
        
    def interpolacaoIterativa(self, x, f, n):
        
        matriz = [[]]
        
        #ordem 0
        l = n
        for i in range (0, l):
            matriz[0].append([i, i, f(x[i])])
           
        #ordem n
        for k in range (0, n-1):
            #valores necessarios e inicializacao da linha da matriz
            l = n-k
            matriz.append([])
            
            for i in range (1, l):
                
                #pega intervalos
                ini = matriz[k][i-1][0]
                fim = matriz[k][i][1]
                
                #posicoes x
                xi = x[ini]
                xf = x[fim]
                
                #posicoes f(x)
                fxi = matriz[k][i-1][2]
                fxf = matriz[k][i][2]
                
                #diferenca dividida
                matriz[k+1].append([ini, fim, self.diferenca(xi, xf, fxi, fxf) ])
        
                if(k == n-2):
                    break;
                    
        coeficientes = []
        for i in range (0, n):
            coeficientes.append(matriz[i][0][2])
            
        print(np.matrix(matriz))
        print(coeficientes)
    

    def interpolacaoRecursiva(self,x, f, n):
        
        #redefine ela como vetor float128
        coeficientes = np.zeros((n,), dtype=np.float128)
        
        #armazena primeiro coeficiente
        coeficientes[0] = x[0]
        
        #comeca a recursao
        self.diferencaRec(x, f, coeficientes)
        
        #coeficientes do polinomio armazenados
        print(coeficientes)
        
    def avaliarPonto(self, x):
        return 0
        
    def diferencaRec(self, x, f, n, c):
          
        valor = 0
        t = len(x)-1

        if(len(x)==2):
            valor = (f(x[1]) - f(x[0])) / (x[1] - x[0])
        else:
            valor = (self.diferencaRec(x[:t], f, n, c) - self.diferencaRec(x[1:t+1], f,n, c)) / (x[t] - x[0])
            
        if(c[t] == 0):
            c[t] = valor
            print(repr("[" + repr(t) + "] = " + repr(valor)))
        
        
        return valor