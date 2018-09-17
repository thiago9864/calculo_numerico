# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 10:24:16 2018

@author: Thiago, Renan
"""
import numpy as np

num_passos_bis = 0
num_passos_fp = 0

### Funcao pra achar a raiz ###

def F(x):
    return x**10 - 1.0
    
## Intervalo ##
    
a = np.dtype('f8')
b = np.dtype('f8')
E = np.dtype('f8')
    
a = 0.0
b = 2.0
ordem = 13
E = 10**-ordem

def insert_dash(string, index):
    return string[:index] + '|' + string[index:]
    
## Métodos ##

def bissecao(_a, _b):

    _x = np.dtype('f8')
    e = _b - _a
    mx = e / 2.0
    _x = mx + _a
    iteracoes = 0
    
    global num_passos_bis

    #checa se a raiz esta nos intervalos
    if(F(_a) == 0.0):
        return _a
    if(F(_b) == 0.0):
        return _b

    while(e > E):

        #acha meio
        e = _b - _a
        mx = e / 2.0
        _x = mx + _a

        #determina valores
        fa = F(_a)
        fb = F(_b)
        fx = F(_x)
        
        #se a raiz estiver entre um dos intervalos formados, ajusta intervalo (_a, _b)
        if((fa >= 0 and fx < 0) or (fa < 0 and fx >= 0)):
            #intervalo a-x
            _b = _x
        else:
            #intervalo x-b
            _a = _x
            
        iteracoes += 1
        
        if(iteracoes > 10000):
            print("saiu no break da bissessao")
            break
    
    num_passos_bis = iteracoes

    #atingiu condicao de parada
    return _x
    
def FalsaPosicao(_a,_b):
    
    x = np.dtype('f8')
    xa = np.dtype('f8')
    
    xa = 0
    x = 0
    iteracoes = 0
    
    global num_passos_fp
    
    #Checa a condição de existencia
    if (F(_a) * F(_b) < 0):
        
        while (abs(F(x)) > E or iteracoes==0):
            
            x = (_a * F(_b) - _b * F(_a)) / (F(_b) - F(_a))
            
            if(F(_a) * F(x) < 0):
                _b = x
            else:
                _a = x
                
            xa = x
            iteracoes += 1
            
            if(iteracoes > 10000):
                print("saiu no break da falsa posicao")
                break
            
        num_passos_fp = iteracoes    
        
        return x
        
        
raiz_bissecao = bissecao(a, b)
raiz_falsa_pos = FalsaPosicao(a, b)

string_bis = repr(raiz_bissecao)
string_fp = repr(raiz_falsa_pos)

print("\nBissecao:\n")
print("Raiz:"+insert_dash(string_bis, ordem + string_bis.find('.') + 1)+"")
print("Precisao:"+str(E)+"")
print("Num de Iteracoes:"+str(num_passos_bis)+"")

print("\nFalsa Posicao:\n")
print("Raiz:"+insert_dash(string_fp, ordem + string_fp.find('.') + 1)+"")
print("Precisao:"+str(E)+"")
print("Num de Iteracoes:"+str(num_passos_fp)+"")