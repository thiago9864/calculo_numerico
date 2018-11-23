# -*- coding: utf-8 -*-
"""
@author: Renan
"""

#Calculo do Erro

import numpy as np
from MetodosIntegracao import NewtonCotes

class Erro():
    
    def erroDaNorma(self,num_part,a,b,f,fi):
      erroL2 = 0    
      '''
      f é a função normal e fi é a função ajustada por minimos quadrados
      '''
      
      for i in range(0,num_part-1):
          # U = f(x) - fi(x)
          erroL2 = erroL2 + NewtonCotes().Simpsom38(,+1,U) 
          '''
          Aqui entram xi e xi+1 para o intervalo de integração e a função a ser integrada
          xi no caso é o intervalo [a,b] e a função f
          '''
      
        erroL2 = np.sqrt(erroL2)
      
      return np.abs(erroL2)

