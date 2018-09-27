# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 10:17:31 2018

@author: Thiago, Renan
"""

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import math as m

### Funcao pra achar a raiz ###

def F(x):
    return m.exp(x) * (x - 1.0) - m.exp(-x) * (x + 1.0)

### Derivada da funcao

def dF(x):
    return m.exp(-x) * (x * m.exp(-2.0 * x) + x)
    
### Função de ponto fixo
    
def pF(x):
    return ((x+1.0)/m.exp(2.0 * x)) + 1.0

  
def insert_dash(string, index):
    return string[:index] + '|' + string[index:]
    

### Intervalo e chute inicial ###
    
a = np.dtype('f8')
b = np.dtype('f8')
xi_s = np.dtype('f8')
xi_n = np.dtype('f8')
E = np.dtype('f8')
    
a = 0.0
b = 2.0
xi = 1.5 #chute inicial
x0 = 1.0 #chute inicial (secante)
ordem = 8
E = 10**-ordem


### Metodos

def bisseccao(_a, _b, Er):

    _xa = _b
    _x = _b
    iteracoes = 0
    lista_iteracoes = []

    while(iteracoes < 10000):

        #acha meio
        e = _b - _a
        mx = e / 2.0
        _x = mx + _a

        #determina valores
        fa = F(_a)
        fx = F(_x)
        
        #se a raiz estiver entre um dos intervalos formados, ajusta intervalo (_a, _b)
        if((fa >= 0 and fx < 0) or (fa < 0 and fx >= 0)):
            #intervalo a-x
            _b = _x
        else:
            #intervalo x-b
            _a = _x
            
        lista_iteracoes.append([iteracoes, _x, e])
        
        if(abs((_x - _xa) / _x) < Er):
            break;
            
        _xa = _x
        iteracoes += 1
        
        
        

    #atingiu condicao de parada
    return [_x, iteracoes, lista_iteracoes]
  
def falsaPosicao(_a,_b, Er):
    
    xa = _b
    x = _b
    iteracoes = 0
    lista_iteracoes = []
    
    #Checa a condição de existencia
    if (F(_a) * F(_b) < 0):
        
        while(iteracoes < 10000):
            
            x = (_a * F(_b) - _b * F(_a)) / (iteracoesF(_b) - F(_a))
            
            if(F(_a) * F(x) < 0):
                _b = x
            else:
                _a = x
                
            
            lista_iteracoes.append([iteracoes, x, abs((x - xa) / x)])
            
            if(abs((x - xa) / x) < Er):
                break;
                
            xa = x            
            iteracoes += 1
        
        
    return [x, iteracoes, lista_iteracoes] 


def pontoFixo(xi, Er):
    _x = xi
    x0 = xi
    iteracoes = 0
    lista_iteracoes = []
    while(iteracoes < 10000):
        _x = pF(_x)
        lista_iteracoes.append([iteracoes, _x, abs((_x - x0) / _x)])
        
        if(abs((_x - x0) / _x) < Er):
            break;
            
        x0 = _x
        iteracoes += 1
        
        
    return [_x, iteracoes, lista_iteracoes]
    
    
def secante(x0, x1, Er):
    iteracoes = 1
    x2 = np.dtype('f8')
    x2 = 0.0
    lista_iteracoes = []
    
    lista_iteracoes.append([0, x1, abs((x1 - x0) / x1) ])
    
    while(iteracoes < 10000):
        f0 = F(x0)
        f1 = F(x1)
        x2 = x2 - ((f1 * (x1 - x0)) / (f1 - f0))
        
        lista_iteracoes.append([iteracoes, x2, abs((x2 - x1) / x2) ])
        
        if(abs((x2 - x1) / x2) < Er):
            break
        
        x0 = x1
        x1 = x2
        iteracoes += 1
        
    return [x2, iteracoes, lista_iteracoes]
    
    
def newton(xi, Er):
    iteracoes = 0
    lista_iteracoes = []
    x0 = xi
    
    while(iteracoes < 10000):
        x1 = x0 - (F(x0) / dF(x0))
        
        lista_iteracoes.append([iteracoes, x1, abs((x1 - x0) / x1)])
        
        if(abs((x1 - x0) / x1) < Er):
            break
        
        x0 = x1
        iteracoes += 1
        
    return [x1, iteracoes, lista_iteracoes]
    
### Implementacao dos metodos ###
    
### Bisseccao
resultado_bisseccao = bisseccao(a, b, E)
raiz_funcao_bis = resultado_bisseccao[0]
num_passos_bis = resultado_bisseccao[1]
string_bis = repr(raiz_funcao_bis)

print("\nBisseccao:\n")
print("Raiz:"+insert_dash(string_bis, ordem + string_bis.find('.') + 1)+"")
print("Precisao:"+str(E)+"")
print("Num de Iteracoes:"+str(num_passos_bis)+"")


### Falsa Posicao
resultado_falsa_pos = falsaPosicao(a, b, E)
raiz_funcao_fpos = resultado_falsa_pos[0]
num_passos_fpos = resultado_falsa_pos[1]
string_fpos = repr(raiz_funcao_fpos)

print("\nFalsa Posicao:\n")
print("Raiz:"+insert_dash(string_fpos, ordem + string_fpos.find('.') + 1)+"")
print("Precisao:"+str(E)+"")
print("Num de Iteracoes:"+str(num_passos_fpos)+"")

### Ponto Fixo
resultado_ponto_fixo = pontoFixo(xi, E)
raiz_funcao_pfixo = resultado_ponto_fixo[0]
num_passos_pfixo = resultado_ponto_fixo[1]
string_pfixo = repr(raiz_funcao_pfixo)

print("\nPonto Fixo:\n")
print("Chute inicial:"+str(xi)+"")
print("Raiz:"+insert_dash(string_pfixo, ordem + string_pfixo.find('.') + 1)+"")
print("Precisao:"+str(E)+"")
print("Num de Iteracoes:"+str(num_passos_pfixo)+"")


### Newton
resultado_newton = pontoFixo(xi, E)
raiz_funcao_newton = resultado_newton[0]
num_passos_newton = resultado_newton[1]
string_newton = repr(raiz_funcao_newton)

print("\nNewton:\n")
print("Chute inicial:"+str(xi)+"")
print("Raiz:"+insert_dash(string_newton, ordem + string_newton.find('.') + 1)+"")
print("Precisao:"+str(E)+"")
print("Num de Iteracoes:"+str(num_passos_newton)+"")

### Secante
resultado_secante = secante(x0, xi, E)
raiz_funcao_secante = resultado_secante[0]
num_passos_secante = resultado_secante[1]
string_secante = repr(raiz_funcao_secante)

print("\nSecante:\n")
print("Chute inicial 0: "+str(x0)+"")
print("Chute inicial 1:"+str(xi)+"")
print("Raiz:"+insert_dash(string_secante, ordem + string_secante.find('.') + 1)+"")
print("Precisao:"+str(E)+"")
print("Num de Iteracoes:"+str(num_passos_secante)+"")


### Geração de gráficos

def geraGrafico(nome, lista):
    lista0 = []
    lista1 = []
    lista2 = []
    
    for item in lista:
        lista0.append(item[0])
        lista1.append(item[1])
        lista2.append(item[2])
        
    #plota grafico
    plt.plot(
        lista0, lista1, 'b--', 
        lista0, lista2, 'r--'
    )
    
    #legendas do grafico
    se_line = mlines.Line2D([], [], color='blue', marker='', markersize=0, label=u'Solução')
    ee_line = mlines.Line2D([], [], color='red', marker='', markersize=0, label=u'Erro')
    
    plt.legend(handles=[se_line, ee_line])
    plt.title(nome)
    plt.show()
    

#geraGrafico(u"Bissecçao", resultado_bisseccao[2])
#geraGrafico(u"Falsa Posicao", resultado_falsa_pos[2])
#geraGrafico(u"Ponto Fixo", resultado_ponto_fixo[2])
#geraGrafico(u"Newton", resultado_newton[2])
geraGrafico(u"Secante", resultado_secante[2])