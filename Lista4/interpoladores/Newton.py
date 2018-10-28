#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 21:09:49 2018

@author: Thiago Almeida
"""
import numpy as np
#from Utils import Utils

class Newton():

    
    '''
    Interpolacao iterativa nao vai ser usada porque so funcionou com lista
    e nao tem tanta precisao quanto a recursiva, por arrays de float128
    '''
    
    '''
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
    ''' 
    
    
    '''
    Solucao do slide, que se comportou da mesma maneira que eu tinha feito
    '''
    
    def interpolacaoSlide(self, n, x, y):
        
        d = np.zeros((n,), dtype=np.float128)
        
        for i in range (0, n): 
            d[i] = y[i]
         
        for k in range (1, n): 
            for i in range (n, k): 
                d[i] = (d[i] - d[i-1]) / (x[i] - x[i-k]);
        
        return d
                
    def avaliarPontoSlide(self, n, z, x, d):
        r = d[n-1]
        
        for i in list(reversed(range (0, n-1))): 
            r = r * (z - x[i]) + d[i]
            
        return r
            
    '''
    Solucao que eu fiz com recursividade
    '''
    
    #constroi o polinomio para interpolacao
    def interpolacao(self,x, y):
        
        global pontos
        global coeficientes
        
        #a ordem depende da quantidade de pontos
        n = len(x)
        
        #armazena pontos
        pontos = np.array(x, dtype=np.float128, copy=True)
        
        #redefine ela como vetor float128
        coeficientes = np.zeros((n,), dtype=np.float128)
        
        #armazena primeiro coeficiente
        coeficientes[0] = y[0]
        
        #comeca a recursao
        self.diferencaRec(x, y, n, coeficientes)
        
        #coeficientes do polinomio armazenados
        #print(list(coeficientes))
        
    #auxiliar da recursao    
    def diferencaRec(self, x, y, n, c):
          
        valor = 0.0
        t = len(x)
        
        #condicao de parada
        if(len(x)==2):
            valor = (y[1] - y[0]) / (x[1] - x[0])
            
        #condicao de recursividade
        elif(len(x)>2):
            tmp_x1 = x[0:t-1]
            tmp_y1 = y[0:t-1]
            tmp_x2 = x[1:t]
            tmp_y2 = y[1:t]
            xi = x[0]
            xf = x[t-1]
            
            valor = (self.diferencaRec(tmp_x1, tmp_y1, n, c) - self.diferencaRec(tmp_x2, tmp_y2, n, c)) / (xf - xi)
        
        #armazena coeficiente
        if(c[t-1] == 0):
            c[t-1] = valor
        
        return valor
    
    #usa os dados da construcao do polinomio pra avaliar os pontos
    def avaliarPonto(self, x):
        
        global pontos
        global coeficientes
        
        n = len(coeficientes)
        valor = 0.0

        for k in range (0, n):
            valor_tmp = coeficientes[k]

            for t in range (0, k):
                valor_tmp *= (x - pontos[t]) 

            valor += valor_tmp

        return valor