"""
@author: Renan
"""

#Calculo do Erro

import numpy as np
from MetodosIntegracao import NewtonCotes

class Erro():
   
    def erroDaNorma(self,K,a,b,f,fi):
      #K = Numero de Partições
      #a = intervalo inicial
      #b = intervalo final
      #f é a função normal e fi é a função ajustada por minimos quadrados
      erroL2 = 0    
        
      
      for i in range(0,K-1):
          U = (f[i]-fi[i])**2
          erroL2 = erroL2 + NewtonCotes().Simpsom38(a,b,K+1,U)
          
      erroL2 = np.sqrt(erroL2) #Calcula a Raiz
      
      return np.abs(erroL2) #Retorna o modulo do resultado
       
         
    
    
    
    
    