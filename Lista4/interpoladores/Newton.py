#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 21:09:49 2018

@author: Thiago Almeida
"""

class Newton(n,x,y):

    def executar(self,n,x,y):
        for i in range(n):
            d[i] = y[i]
        for k in range (1,n):
            d[i] = (d[i] - d[i-1])/(x[i] - x[i-k])
        return d