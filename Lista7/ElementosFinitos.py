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
    
    def funcaoBase(self):
        