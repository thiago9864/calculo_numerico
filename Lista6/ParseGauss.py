#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 25 19:56:52 2018

@author: thiagoalmeida
"""
#import numpy as np

inicio_ordem = 2
fim_ordem = 64
limite_casas_decimais = 16


##### Parse da matriz dos t

arquivo = open("lgvalues-abscissa.txt", 'r');

conteudo = arquivo.readlines();
string = ""
for i in conteudo:
   string += i

arquivo.close() 

novo_arquivo = "t = [\n"

for k in range(inicio_ordem, fim_ordem+1):

    num_linhas = k
    stri = '$legendre_roots['+str(num_linhas)+'] = array('
    strf = ');'
    ind_i = string.find(stri) + len(stri)
    ind_f = string.find(strf, ind_i)
    dados = string[ind_i:ind_f].replace('\t', '').replace('\r', '').replace('\n', '').replace(' ', '')
    dados_arr_s = dados.split(',')
    
    novo_arquivo += '['
    
    
    for i in range(num_linhas):
        if(i != 0):
            novo_arquivo += ','
            
        numero = dados_arr_s[i]
        ponto = numero.find('.')
        corte = ponto + limite_casas_decimais + 1
        novo_arquivo += dados_arr_s[i][:corte]
    
    if(k != fim_ordem):
        novo_arquivo += '],\n'
    else:
        novo_arquivo += ']\n'

novo_arquivo += "]\n"
        
        
##### Parse da matriz dos w
        
        
        
arquivo2 = open("lgvalues-weights.txt", 'r');

conteudo = arquivo2.readlines();
string = ""
for i in conteudo:
   string += i

arquivo2.close() 

novo_arquivo += "w = [\n"

for k in range(inicio_ordem, fim_ordem+1):

    num_linhas = k
    stri = '$quadrature_weights['+str(num_linhas)+'] = array('
    strf = ');'
    ind_i = string.find(stri) + len(stri)
    ind_f = string.find(strf, ind_i)
    dados = string[ind_i:ind_f].replace('\t', '').replace('\r', '').replace('\n', '').replace(' ', '')
    dados_arr_s = dados.split(',')
    
    novo_arquivo += '['
    
    
    for i in range(num_linhas):
        if(i != 0):
            novo_arquivo += ','
            
        numero = dados_arr_s[i]
        ponto = numero.find('.')
        corte = ponto + limite_casas_decimais + 1
        novo_arquivo += dados_arr_s[i][:corte]
    
    if(k != fim_ordem):
        novo_arquivo += '],\n'
    else:
        novo_arquivo += ']\n'

novo_arquivo += "]"
        

print(novo_arquivo)

file = open("testfile.py","w") 
 
file.write(novo_arquivo) 
 
file.close() 
    