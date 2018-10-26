#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 21:05:05 2018

@author: Renan Nunes
"""

class Lagrange():

    def executar(self,n,x,y,z):
        
        r = 0
        for i in range (1,n):
            c = 1
            d = 1
            for j in range (1,n):
                if (i != j):
                    c = c * (z - x[j])
                    d = d * (x[i] - x[j])
            r = r + y[i] * (c/d)
        
        return r