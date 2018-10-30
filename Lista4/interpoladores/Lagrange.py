#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 21:05:05 2018

@author: Renan Nunes
"""

class Lagrange():

    #constroi o polinomio para interpolacao
    def interpolacao(self,n,x,y,z):
        
        r = 0
        for i in range (0,n):
            c = 1.0
            d = 1.0
            for j in range (0,n):
                if (i != j):
                    c = c * (z - x[j])
                    d = d * (x[i] - x[j])
            r = r + y[i] * (c/d)
        
        return r
    