# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import matplotlib.lines as mlines

arquivo = open("TesteData.txt", 'r');
datalist = arquivo.readlines() 

tempo = [0] * len(datalist)
solucao_exata = [0] * len(datalist)
solucao_aproximada = [0] * len(datalist)
metodo = "C++ Gauss"

listsize = len(datalist)

for i in range(listsize):
    data = datalist[i].split("|")
    #print(data)
    tempo[i] = float(data[0])
    solucao_exata[i] = float(data[1])
    solucao_aproximada[i] = float(data[2].replace("\n", ""))

print("len tempo: ", len(tempo))
print("len solucao_aproximada: ", len(solucao_aproximada))
print("len solucao_exata: ", len(solucao_exata))

#insere intervalos de contorno
#solucao_aproximada.insert(0, 0)
#solucao_aproximada.append(0)

 #plota grafico da função
plt.plot(
    tempo, solucao_exata, 'b--',
    tempo, solucao_aproximada, 'r-'    
    )
plt.ylabel(u"Valor de u(h)") #esse 'u' antes da string é pra converter o texto pra unicode
plt.xlabel(u"Valor de h, " + str(len(tempo)-1) + u" partições")

#legendas do grafico
se_line = mlines.Line2D([], [], color='blue', marker='', markersize=0, label=u'Sol Exata')
ac_line = mlines.Line2D([], [], color='red', marker='', markersize=0, label=u'Met. ' + metodo)

plt.legend(handles=[se_line, ac_line])

plt.title(u"Metodo de "+metodo+u" x Solução exata", )

#plt.axis([0, 50, 0, 100])
plt.show() 