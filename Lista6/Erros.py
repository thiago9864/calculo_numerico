# -*- coding: utf-8 -*-

"""
@author: Renan, Thiago
"""

#Calculo do Erro

import numpy as np
from MetodosIntegracao import NewtonCotes
from MinimosQuadrados import MinimosQuadrados

class Erro():
   
   
    def erroDaNorma(self,K,a,b,f,coeficientes, n, particoes): #1
      #K = Numero de Particoes
      #a = intervalo inicial
      #b = intervalo final
      #f = funcao normal
      #coeficientes, n, particoes = necessario pra funcao fi
    
      erroL2 = 0    
      
      #funcao pra jogar na integral
      def U(x):
          fi = MinimosQuadrados().interpolaCoeficientes(coeficientes, n, x)
          return (f(x)-fi)**2;
      
      for i in range(0,K-1):
          #intervalo de integracao
          xi = particoes[i]
          xf = particoes[i+1]
          
          erroL2 = erroL2 + NewtonCotes().Simpsom38(xi,xf,U)#integracao numerica
          
      erroL2 = np.sqrt(erroL2) #Calcula a Raiz
      
      return np.abs(erroL2) #Retorna o modulo do resultado
       
         
    
    
    
    
    