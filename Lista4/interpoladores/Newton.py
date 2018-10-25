#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 21:09:49 2018

@author: Thiago Almeida
"""
import numpy as np
import math as m

class Newton():
    
    def executar(self,x, f, n):
        return self.diferenca(x, f)
        
    def diferenca(self, x, f):
        if(len(x)==2):
            return (f(x[1]) - f(x[0])) / (x[1] - x[0])
        else:
            t = len(x)-1
            return (self.diferenca(x[:t], f) - self.diferenca(x[1:t+1], f)) / (x[t] - x[0])
    