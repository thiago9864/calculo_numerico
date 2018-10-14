# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import matplotlib.lines as mlines

#imports locais
from Utils import Utils
from metodos_numericos.Gauss import Gauss
from metodos_numericos.LU import LU
    
#### Gera matriz para essa questao #####
    
    
def gerarMatriz(ordem):
    
    M = Utils().inicializaMatriz(ordem)
    B = [0] * ordem
    
    #cria a matriz
    for i in range(ordem):
        for j in range(ordem):
            M[i][j] = 1.0 / (i + j + 1.0)
    
    #cria o vetor fonte
    for i in range(ordem):
        B[i] = 1.0 / (i + ordem + 1.0)
        
    return [M, B]
   
##### Gerador de grafico #####

def gerarGrafico(tempo, solucao_aproximada, titulo, metodo):
    
    #print("len tempo: ", len(tempo))
    #print("len solucao_aproximada: ", len(solucao_aproximada))
    #print("len solucao_exata: ", len(solucao_exata))
    
    
    #plota grafico da função
    plt.plot(
        tempo, solucao_aproximada, 'r-'    
        )
    plt.ylabel(u"Valor de u(h)") #esse 'u' antes da string é pra converter o texto pra unicode
    plt.xlabel(u"Valor de h, " + str(len(tempo)) + u" elementos")
    
    #legendas do grafico
    ac_line = mlines.Line2D([], [], color='red', marker='', markersize=0, label=u'Met. ' + metodo)
    
    plt.legend(handles=[ac_line])
    
    plt.title(titulo)
    
    #plt.axis([0, 50, 0, 100])
    plt.show()   
     
#### Executa os metodos #####

numero_de_elementos = 1000
    
mat = gerarMatriz(numero_de_elementos)

M = mat[0]
B = mat[1]

#gera eixo-X do grafico
arrEixoX = []
for i in range(numero_de_elementos):
    arrEixoX.append(i)

#no mac executa aproximadamente 3837155 operacoes por segundo
#formula para o tempo previsto no mac -> ((numero_de_elementos^3)*1/3) / 3837155 = segundos_rodando
#Utils().imprimeMatriz(M, B)
#Utils().obtemInfoMatriz(M)


print("Metodo de Gauss (sem pivoteamento)")
inicio = Utils().getTime()
resGauss = Gauss().executar(M, B)
fim  = Utils().getTime()
print("Tamanho da matriz: " + repr(numero_de_elementos) + "x" + repr(numero_de_elementos))
print("Passos ate a resolucao: " + repr(resGauss[1]))
Utils().imprimeDiferencaTempo(inicio, fim)
#print("Resultado do sistema por Gauss: ", resGauss[0]);
#gerarGrafico(arrEixoX, resGauss[0], "Metodo de Gauss", "Gauss")
erro = Utils().erroResidual(M, resGauss[0], B)
print("Erro maximo encontrado: " + repr(erro[1]))
gerarGrafico(arrEixoX, erro[0], "Erro no Metodo de Gauss", "Gauss")

print("")
print("------")
print("")

print("Metodo de Gauss (com pivoteamento)")
inicio = Utils().getTime()
resGaussP = Gauss().executarComPivoteamento(M, B)
fim  = Utils().getTime()
print("Tamanho da matriz: " + repr(numero_de_elementos) + "x" + repr(numero_de_elementos))
print("Passos ate a resolucao: " + repr(resGaussP[1]))
Utils().imprimeDiferencaTempo(inicio, fim)
#print("Resultado do sistema por Gauss Pivoteado: ", resGaussP[0]);
#gerarGrafico(arrEixoX, resGaussP[0], "Metodo de Gauss", "Gauss (com pivoteamento)")
erro = Utils().erroResidual(M, resGaussP[0], B)
print("Erro maximo encontrado: " + repr(erro[1]))
gerarGrafico(arrEixoX, erro[0], "Erro no Metodo de Gauss (Pivoteado)", "Gauss")

print("")
print("------")
print("")

print("Metodo LU")
inicio = Utils().getTime()
resLU = LU().executar(M, B)
fim  = Utils().getTime()
print("Tamanho da matriz: " + repr(numero_de_elementos) + "x" + repr(numero_de_elementos))
print("Passos ate a resolucao: " + repr(resLU[1]))
Utils().imprimeDiferencaTempo(inicio, fim)
#print("Resultado do sistema por LU: ", resLU[0]);
#gerarGrafico(arrEixoX, resLU[0], "Metodo LU", "LU")
erro = Utils().erroResidual(M, resLU[0], B)
print("Erro maximo encontrado: " + repr(erro[1]))
gerarGrafico(arrEixoX, erro[0], "Erro no Metodo LU", "LU")