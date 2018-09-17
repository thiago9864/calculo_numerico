# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 11:40:34 2018

@author: Thiago, Renan
"""
import numpy as np
import math as m

num_passos_f3 = 0
num_passos_f4 = 0

def insert_dash(string, index):
    return string[:index] + '|' + string[index:]
    

### Funcao pra achar a raiz ###

def F(x):
    return x**3 + 4.0 * x**2 - 10.0
    
### Funcoes de ponto fixo que convergem, de acordo com a demonstracao ###
    
def f3(x):
    return m.sqrt(10.0 / (4.0+x))
    
def f4(x):
    return x - ((x**3 + 4.0*x**2 - 10.0) / (3.0 * x**2 + 8.0*x))
    
### Intervalo e chute inicial ###
    
a = np.dtype('f8')
b = np.dtype('f8')
xi = np.dtype('f8')
E = np.dtype('f8')
    
a = 1.0
b = 2.0
xi = 1.5 #chute inicial

ordem = 5
E = 10**-ordem


### Implementação do metodo

def pontoFixo(xi, funcao):
    _x = xi
    iteracoes = 0
    while(abs(F(_x)) > E):
        
        _x = funcao(_x)
        iteracoes += 1
        
        if(iteracoes >= 10000):
            print("saiu no break")
            break
        
    return [_x, iteracoes]
    
### Obtem os resultados
    
p3 = pontoFixo(xi, f3)
p4 = pontoFixo(xi, f4)

raiz_funcao3 = p3[0]
num_passos_f3 = p3[1]
raiz_funcao4 = p4[0]
num_passos_f4 = p4[1]

string_3 = repr(raiz_funcao3)
string_4 = repr(raiz_funcao4)

print("\nFuncao 3:\n")
print("Raiz:"+insert_dash(string_3, ordem + string_3.find('.') + 1)+"")
print("F(x):"+repr(F(raiz_funcao3)))
print("Precisao:"+str(E)+"")
print("Num de Iteracoes:"+str(num_passos_f3)+"")

print("\nFuncao 4:\n")
print("Raiz:"+insert_dash(string_4, ordem + string_4.find('.') + 1)+"")
print("F(x):"+repr(F(raiz_funcao4)))
print("Precisao:"+str(E)+"")
print("Num de Iteracoes:"+str(num_passos_f4)+"")