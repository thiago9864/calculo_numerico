# -*- coding: utf-8 -*- 
"""
Created on Mon Sep 17 16:00:20 2018

@author: Thiago, Renan
"""
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import math as m

### Funcao pra achar a raiz ###

def F(x):
    return m.exp(2.0*x) - 2.0 * m.exp(x) + 1

### Derivada da funcao

def dF(x):
    return 2 * (m.exp(x) - 1) * m.exp(x)

  
def insert_dash(string, index):
    return string[:index] + '|' + string[index:]
   

### Intervalo e chute inicial ###
    
a = np.dtype('f8')
b = np.dtype('f8')
xi_s = np.dtype('f8')
xi_n = np.dtype('f8')
E = np.dtype('f8')
    
a = -1.0
b = 1.0
xi_s = 0 #chute inicial (secante)
xi_n = 1 #chute inicial (Newton)

ordem = 5
E = 10**-ordem


### Metodos

def bissecao(_a, _b, Er):

    _x = np.dtype('f8')
    e = _b - _a
    mx = e / 2.0
    _x = mx + _a
    iteracoes = 0

    #checa se a raiz esta nos intervalos
    if(F(_a) == 0.0):
        return _a
    if(F(_b) == 0.0):
        return _b

    while(e > Er):

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
            
        iteracoes += 1
        
        if(iteracoes > 10000):
            print("saiu no break da bissessao")
            break

    #atingiu condicao de parada
    return [_x, iteracoes]
    
    
def secante(x0, x1):
    iteracoes = 0
    x2 = np.dtype('f8')
    x2 = 0.0
    lista_iteracoes = []
    
    while(iteracoes < 10000):
        f0 = F(x0)
        f1 = F(x1)
        x2 = x2 - ((f1 * (x1 - x0)) / (f1 - f0))
        
        lista_iteracoes.append([iteracoes, x2, abs(x2 - x1)])
        
        if(abs(x2 - x1) < E):
            break
        
        x0 = x1
        x1 = x2
        iteracoes += 1
        
    return [x2, iteracoes, lista_iteracoes]
    
    
def newton(xi):
    iteracoes = 0
    lista_iteracoes = []
    x0 = xi
    
    while(iteracoes < 10000):
        x1 = x0 - (F(x0) / dF(x0))
        
        lista_iteracoes.append([iteracoes, x1, abs(x1 - x0)])
        
        if(abs(x1 - x0) < E):
            break
        
        x0 = x1
        iteracoes += 1
        
    return [x1, iteracoes, lista_iteracoes]
    
    
### Executa metodo da secante com primeiro chute
    
lista0 = []
lista1 = []
lista2 = []

#obtem primeiro chute
xi_b = bissecao(a, b, 0.1)

resultado_secante = secante(xi_s, xi_b[0])

raiz_funcao_sec = resultado_secante[0]
num_passos_sec = resultado_secante[1]
lista_iteracoes_sec = resultado_secante[2]

string_sec = repr(raiz_funcao_sec)

print("\nFuncao com chute inicial "+repr(xi_s)+" e metodo da Secante:\n")
print("Raiz:"+insert_dash(string_sec, ordem + string_sec.find('.') + 1)+"")
print("F(x):"+repr(F(raiz_funcao_sec)))
print("Precisao:"+str(E)+"")

print("Iteracoes")
print("#".ljust(5) + "|" + "valor de x".ljust(25) + "|" + "erro".ljust(25) )
for item in lista_iteracoes_sec:
    lista0.append(item[0])
    lista1.append(item[1])
    lista2.append(item[2])
    print(repr(item[0]).ljust(5) + "|" + repr(item[1]).ljust(25) + "|" + repr(item[2]).ljust(25) )
    
print("Num de Iteracoes:"+str(num_passos_sec)+"")

#plota grafico
plt.plot(
    lista0, lista1, 'b--', 
    lista0, lista2, 'r--'
)

#legendas do grafico
se_line = mlines.Line2D([], [], color='blue', marker='', markersize=0, label=u'Solução')
ee_line = mlines.Line2D([], [], color='red', marker='', markersize=0, label=u'Erro')

plt.legend(handles=[se_line, ee_line])
plt.title("Metodo da Secante")
plt.show()

### Executa metodo de Newton com o segundo chute

lista0 = []
lista1 = []
lista2 = []

resultado_newton = newton(xi_n)

raiz_funcao_new = resultado_newton[0]
num_passos_new = resultado_newton[1]
lista_iteracoes_new = resultado_newton[2]

string_new = repr(raiz_funcao_new)

print("\nFuncao com chute inicial "+repr(xi_n)+" e metodo de Newton:\n")
print("Raiz:"+insert_dash(string_new, ordem + string_new.find('.') + 1)+"")
print("F(x):"+repr(F(raiz_funcao_new)))
print("Precisao:"+str(E)+"")

print("Iteracoes")
print("#".ljust(5) + "|" + "valor de x".ljust(25) + "|" + "erro".ljust(25) )
for item in lista_iteracoes_new:
    lista0.append(item[0])
    lista1.append(item[1])
    lista2.append(item[2])
    print(repr(item[0]).ljust(5) + "|" + repr(item[1]).ljust(25) + "|" + repr(item[2]).ljust(25) )
    
print("Num de Iteracoes:"+str(num_passos_new)+"")

#plota grafico
'''
plt.plot(
    lista0, lista1, 'b--', 
    lista0, lista2, 'r--'
)

#legendas do grafico
se_line = mlines.Line2D([], [], color='blue', marker='', markersize=0, label=u'Solução')
ee_line = mlines.Line2D([], [], color='red', marker='', markersize=0, label=u'Erro')

plt.legend(handles=[se_line, ee_line])
plt.title("Metodo de Newton")
plt.show()
'''