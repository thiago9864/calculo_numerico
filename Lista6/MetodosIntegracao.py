# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 16:38:23 2018

@author: thiagoalmeida
"""

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
        
    def Retangulo(self):
        return 0
        
    def PontoMedio(self):
        return 0
        
    def Trapezio(self):
        return 0
        
    def Simpsom13(self, x, f):
        tam = len(x)
        
        a = x[0]
        b = x[tam-1]
        h = (b - a) / 2.0
        
        soma = 0
        for i in range (0, len(x)):
            
            if(i==0 or i==tam-1):
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
        
    def Simpsom38(self, x, f):
        tam = len(x)
        
        a = x[0]
        b = x[tam-1]
        h = (b - a) / 3.0
        
        soma = 0
        for i in range (0, len(x)):
            
            if(i==0 or i==tam-1):
                #extremos, igualar a 1
                ci = 1.0
            elif(i % 4 == 0):
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