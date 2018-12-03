# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 11:22:02 2018

@author: thiagoalmeida
"""

from scipy.interpolate import lagrange
import numpy as np
from MetodosIntegracao import NewtonCotes
#from MinimosQuadrados import MinimosQuadrados

class ElementosFinitos():
    
    '''
    Fazer uma tabela de funcoes base e derivadas delas
    olhar -> https://jayemmcee.wordpress.com/lagrange-polynomial-interpolation/
    '''
        
    
    def funcaoBase(self, n, x, xk):
        
        L = np.zeros((n,), dtype=np.float128)
        dL = np.zeros((n,), dtype=np.float128)
        
        #funcao de base
        for i in range(n):
            c = 1.0
            d = 1.0
            for j in range(n):
                if(i != j):
                    c = c * (xk - x[j])
                    d = d * (x[i] - x[j])
            L[i] = c / d
            
        #derivada da funcao de base
        for i in range(n):
            c = 0
            d = 1.0
            print("i=" + str(i))
            
            for j in range(n):
                if(i != j):  
                    d = d * (x[i] - x[j])
            
            #calcula as derivadas pro regra do produto  
            for j in range(n):
                if(i != j):
                    print("j=" + str(j))
                    #calcula regra do produto
                    r = 1.0
                    for d in range(n):
                       if(d != j and d != i):
                           print("d=" + str(d))
                           r = r * (xk - x[d])
                       #soma regra do produto
                       c += r
                   
            print("-")             
            dL[i] = c / d
            
        return (L, dL)
                    
            
            