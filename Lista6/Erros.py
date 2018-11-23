# -*- coding: utf-8 -*-
"""
@author: Renan
"""

#Calculo do Erro

import numpy as np
from MetodosIntegracao import NewtonCotes

class Erro():
    
    def erroDaNorma(self,num_part):
      erroL2 = 0    
      
      for i in range(0:num_part-1):
          erroL2 = erroL2 + NewtonCotes().Simpsom38(,+1,)
      
        erroL2 = np.sqrt(erroL2)
      
      return np.abs(erroL2)

