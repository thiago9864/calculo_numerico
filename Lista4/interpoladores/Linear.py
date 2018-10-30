#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 19:24:25 2018

@author: Thiago Almeida
"""

class Linear():

    #constroi o polinomio para interpolacao
    def interpolacao(self,n,x,y,z):
        
        #escolhe dois pontos onde a < z <b
        a=0
        b=0
        for i in range (1, n):
            a = i-1
            b = i
            
            if(x[a] <= z and z <= x[b]):
                break
            
        y0 = y[a]
        x0 = x[a]
        y1 = y[b]        
        x1 = x[b]
        
        return y0 + ((y1 - y0) / (x1 - x0)) * (z - x0)