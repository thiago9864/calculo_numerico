# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 11:22:02 2018

@author: thiagoalmeida
"""

import numpy as np
from TabelaGaussLegendre import TabelaGaussLegendre
from metodos_numericos.LU import LU
from metodos_numericos.Gauss import Gauss
from metodos_numericos.Thomas import Thomas


class ElementosFinitos:
    
    '''
    construtor
    a, b: intervalo [a, b]
    c1, c2: condicoes de contorno
    elementos: numero de elementos a avaliar
    grau_polinomio: grau do polinomio interpolador (pode ir de 1 a 63)
    E: o epsilon
    dataType: tipo de dado. vai ser np.float128 ou np.float64
    '''
    def __init__(self, a, b, y1, y2, elementos, grau_polinomio, E, dataType):
        
        #parametros recebidos
        self.a = a
        self.b = b
        self.y1 = y1
        self.y2 = y2
        self.elementos = elementos
        self.grau_polinomio = grau_polinomio
        self.E = E
        self.dataType = dataType

        
        #parametros calculados
        self.pontos_elementos = int(elementos + 1)
        self.pontos_polinomio = int(grau_polinomio + 1)
        self.tamanho_problema = int(grau_polinomio * elementos + 1)
        
        #inicia vetor dos pontos do intervalo de elementos
        self.x = np.zeros((self.pontos_elementos,), dtype=self.dataType)
        
        #calcula largura do elemento
        self.h = (b-a)/elementos
        
        #gera vetor com pontos do intervalo de elementos
        d = a
        for i in range(self.pontos_elementos):
            self.x[i] = d
            d += self.h
        
        #debug
        #print("h:", self.h)
        #print("x:", self.x)
    
    '''
    Funcao base
    i = indice da funcao fi desejada
    n = numero de pontos    
    xv = vetor com os pontos de x do intervalo 
    xk = ponto a avaliar
    '''
    def fi(self, i, xv, xk):
        
        L = 0.0
        dL = 0.0
            
        #armazena os indices do vetor x pra usar na derivada
        ind = []
        
        '''
        print("--- fi ---")
        print("i", i)  
        print("xv", xv)  
        print("xk", xk)  
        '''
        
        #calcula o valor da funcao de base
        c = 1.0
        d = 1.0
        for j in range(self.pontos_polinomio):
            if(i != j):
                ind.append(j)
                c = c * (xk - xv[j])
                d = d * (xv[i] - xv[j])
        L = c / d
        
        #calcula o valor da derivada da funcao de base
        soma = 0.0
        for j in range(self.pontos_polinomio-1):
            k = ind[j]
            
            if(self.pontos_polinomio > 2):
                soma = soma + (xk - xv[k]) / d
            else:
                soma = soma + 1.0 / d
        
        dL = soma   
        
        '''
        print("d", d)
        print("L", L)  
        print("dL", dL)  
        '''
        
        return (L, dL)
    
   
    '''
    calcula o valor de K da matriz
    i = linha
    j = coluna
    h = tamanho do elemento
    pontos = numero de pontos do polinomio
    '''
    def K(self, i, j):
        
        #print("i: " + str(i) + ", j: " + str(j))
        
        #recupera pontos da tabela
        tw = TabelaGaussLegendre().getValores(self.pontos_polinomio)
        
        #vetor local dos pontos de integracao
        xp = np.zeros((self.pontos_polinomio,), dtype=self.dataType)
        
        a = -1.0
        b = 1.0
        hp = (b-a) / self.grau_polinomio
        d = a
        
        #gera pontos de integração
        for k in range(self.pontos_polinomio):
            xp[k] = d
            d += hp
            
        
        #calcula a integral por Gauss
        soma = 0.0
        for k in range (0, self.pontos_polinomio):   
            ti = self.dataType(tw[0][k])
            wi = self.dataType(tw[1][k])
            
            xk = ti
            
            (f1, df1) = self.fi(i, xp, xk)
            (f2, df2) = self.fi(j, xp, xk)

            #slide
            #f = df1 * (2.0/self.h) * df2 * (2.0/self.h)
            
            #trabalho
            f = df1 * (2.0/self.h) * df2 * (2.0/self.h) * self.E + (f1 * f2)
            
            soma += wi * f * (self.h/2.0)
            
            '''
            #debug
            print("xp", xp)
            print("ti", ti)
            print("wi", wi)
            print("df1", df1)
            print("df2", df2)
            print("h", self.h)
            print("-------")
            '''
        return soma
        
    
    '''
    calcula o valor de F do vetor fonte
    i = linha
    h = tamanho do elemento
    pontos = numero de pontos do polinomio
    '''
    def F(self, i):
        
        #print("i: " + str(i))
        
        #recupera pontos da tabela
        tw = TabelaGaussLegendre().getValores(self.pontos_polinomio)
        
        #gera pontos
        xp = np.zeros((self.pontos_polinomio,), dtype=self.dataType)
        
        q = 2.0 / self.grau_polinomio
        d = -1.0
        
        for k in range(self.pontos_polinomio):
            xp[k] = d
            d += q
            
        
        #calcula a integral por Gauss
        soma = 0.0
        for k in range (0, self.pontos_polinomio):   
            ti = self.dataType(tw[0][k])
            wi = self.dataType(tw[1][k])
            
            xk = ti
            
            (f1, df1) = self.fi(i, xp, xk)
            
            f = f1 * 1.0 #esse 1 e a funcao f do slide
            
            soma += wi * f * (self.h/2.0)
            
            #print("xp", xp)
            #print("ti", ti)
            #print("wi", wi)
            #print("f1", f1)
            #print("h", h)
            #print("-------")
            
        return soma
        
        
    '''
    gera matriz local
    '''
    def matrizLocal(self):
        K = np.zeros((self.pontos_polinomio, self.pontos_polinomio), dtype=self.dataType)
        for i in range(self.pontos_polinomio):
            for j in range(self.pontos_polinomio):
                K[i][j] = self.K(i, j)
        return K
        
    '''
    gera vetor local
    '''
    def vetorLocal(self):
        F = np.zeros((self.pontos_polinomio,), dtype=self.dataType)
        for i in range(self.pontos_polinomio):
                F[i] = self.F(i)
        return F
     
    
    '''
    metodo que preenche uma matriz maior com uma menor a partir de um par de coordenadas
    M: matriz maior
    A: matriz menor
    i, j: coordenadas iniciais
    '''
    def preencheMatrizEmMatriz(self, M, A, i, j):
        ordem_da_menor = len(A[0])
        for im in range(ordem_da_menor):
            for jm in range(ordem_da_menor):
                if(im == 0 and jm == 0):
                    #se for o indice (0, 0), soma com o valor que ja esta la
                    M[im+i][jm+j] += A[im][jm]
                else:
                    #se for qualquer outro, ja substitui direto
                    M[im+i][jm+j] = A[im][jm]
    '''
    metodo que preenche um vetor maior com um vetor menor a partir de um indice
    v: vetor maior
    a: vetor menor
    i: indice de partida
    '''        
    def preencheVetorEmVetor(self, v, a, i):
        ordem_vetor_menor = len(a)
        for im in range(ordem_vetor_menor):
            if(im == 0):
                #se o elemento estiver na area de sobreposicao, soma com o valor que ja esta la
                v[im+i] += a[im]
            else:
                #se for qualquer outro, ja substitui direto
                v[im+i] = a[im]
                
    '''
    gera matriz final
    '''       
    def matrizRigidez(self):
        Kr = np.zeros((self.tamanho_problema, self.tamanho_problema), dtype=self.dataType)
        K = self.matrizLocal()
        
        #preenche matriz de rigidez com as matrizes locais
        for m in range(self.elementos):
            p = m * self.grau_polinomio
            self.preencheMatrizEmMatriz(Kr, K, p, p)
        
        
        #zera contorno da matriz
        for i in range(self.tamanho_problema):
            for j in range(self.tamanho_problema):
                if(i == 0 or i == (self.tamanho_problema-1) or j == 0 or j == (self.tamanho_problema-1)):
                    Kr[i][j] = 0
        
        #impoe condicoes de contorno          
        Kr[0][0] = 1.0
        Kr[self.tamanho_problema-1][self.tamanho_problema-1] = 1.0

                
        return Kr
    
    '''
    gera vetor final
    '''
    def vetorForca(self):
        Fr = np.zeros((self.tamanho_problema,), dtype=self.dataType)
        F = self.vetorLocal()
        
        #preenche vetor forca com os vetores locais
        for m in range(self.elementos):
            p = m * self.grau_polinomio
            self.preencheVetorEmVetor(Fr, F, p)
        
        #impoe condicoes de contorno
        Fr[0] = 0.0
        Fr[self.tamanho_problema-1] = 0.0
        
        return Fr
                
    '''
    metodo principal que calcula os pontos
    '''           
    def calcular(self):
        Kr = self.matrizRigidez()
        Fr = self.vetorForca()
        
        if(self.grau_polinomio == 1):
            print("Resolvendo com algoritmo de Thomas")
            u = Thomas().executar(Kr, Fr)
        else:
            print("Resolvendo com algoritmo LU")
            u = LU().executar(Kr, Fr, self.dataType)
            
            #print("Resolvendo com algoritmo Gauss com pivoteamento")
            #u = Gauss().executarComPivoteamento(Kr, Fr, self.dataType)
        
        return u
                
                