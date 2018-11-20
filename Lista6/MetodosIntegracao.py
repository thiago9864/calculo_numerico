# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 16:38:23 2018

@author: thiagoalmeida
"""
import numpy as np

class NewtonCotes():
        
    def Retangulo(self, a, b, f):
        h = b - a
        x0 = b
        return h * f(x0)
        
    def PontoMedio(self, a, b, f):
        x0 = (b - a) / 2.0
        h = b - a
        return h * f(x0)
        
    def Trapezio(self, a, b, f):
        h = (b - a)
        x0 = a
        x1 = b
        return (h / 2.0) * (f(x0)+f(x1))      
        
    def Simpsom13(self, a, b, f):
        h = (b - a) / 2.0
        x0 = a
        x1 = h
        x2 = b
        return (h / 3.0) * (f(x0) + (4.0 * f(x1)) + f(x2))
        
    def Simpsom38(self, a, b, f):
        h = (b - a) / 3.0
        x0 = a
        x1 = x0 + h
        x2 = x1 + h
        x3 = b
        return ((3.0 * h) / 8.0) * (f(x0) + (3.0 * f(x1)) + (3.0 * f(x2)) + f(x3))
            
            
class Repetidos():
    
    def montaVetorX(self, a, b, m):      
        dm = (b-a)/(m-1)
        x = np.zeros((m,), dtype=np.float128)
        for i in range (0, m):
            x[i] = dm * i
        return x
        
    def Retangulo(self, a, b, m, f):
        x = self.montaVetorX(a, b, m)
            
        h = (b - a) / m
        
        soma = 0
        for i in range (1, len(x)):
            soma += h * f(x[i-1])
            
        return soma
        
    def PontoMedio(self, a, b, m, f):
        x = self.montaVetorX(a, b, m)
            
        h = (b - a) / m
        
        soma = 0
        for i in range (1, len(x)):
            soma += h * f((x[i-1] + x[i-1]) / 2.0)
            
        return soma
        
    def Trapezio(self, a, b, m, f):

        x = self.montaVetorX(a, b, m)
            
        h = (b - a) / m
        
        soma = 0
        for i in range (0, len(x)):
            
            if(i==0 or i==m-1):
                #extremos, igualar a 1
                ci = 1.0
            else:
                #exceto os extremos
                ci = 2.0

                
            soma += f(x[i]) * ci
            
        return soma * (h/2.0) 
        
    def Simpsom13(self, a, b, m, f):
        
        x = self.montaVetorX(a, b, m)
            
        h = (b - a) / m
        
        soma = 0
        for i in range (0, len(x)):
            
            if(i==0 or i==m-1):
                #extremos, igualar a 1
                ci = 1.0
            elif(i % 2 == 0):
                #se for par
                ci = 2.0
            else:
                #se for impar
                ci = 4.0
                
            soma += f(x[i]) * ci
            
        return soma * (h/3.0) 
        
    def Simpsom38(self, a, b, m, f):

        x = self.montaVetorX(a, b, m)
            
        #define h 
        h = (b - a) / m
        
        soma = 0
        for i in range (0, len(x)):
            
            if(i==0 or i==m-1):
                #extremos, igualar a 1
                ci = 1.0
            elif(i % 3 == 0):
                #de quatro em quatro termos e 2
                ci = 2.0
            else:
                #nos demais termos e 3
                ci = 3.0
                
            soma += f(x[i]) * ci
            
        return soma * ((3.0*h)/8.0)            
            
            
class QuadraturaGauss():
    
    def teste(self):
        return 0