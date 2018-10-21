#include <iostream>
#include <cmath>        // std::abs
#include <fstream> 
#include <string>

using namespace std;

#define NUM_ELEMENTOS 999
/*
void imprimeMatriz (double m[NUM_ELEMENTOS][NUM_ELEMENTOS])
{
    int i, j;
    for(i = 0; i < 6; i++)
    {
       for(j = 0; j < 7; j++)
       {
          printf("%10f", m[i][j]);
       }
       printf("\n");
    }
}
*/
void trocarLinhas(double M[NUM_ELEMENTOS][NUM_ELEMENTOS], double B[NUM_ELEMENTOS], int linha_a, int linha_b){
    
	int j;
	int tam = NUM_ELEMENTOS;
	double elemento;
	for (j=0; j<tam; j++){
		//troca a linha da matriz
		elemento = M[linha_a][j];
		M[linha_a][j] = M[linha_b][j];
		M[linha_b][j] = elemento;
	}

	//troca a linha do vetor fonte
	elemento = B[linha_a];
	B[linha_a] = B[linha_b];
	B[linha_b] = elemento;
}

void gerarMatrizQuestao2(double M[NUM_ELEMENTOS][NUM_ELEMENTOS], double B[NUM_ELEMENTOS]){
    
	int i, j;
	int tam = NUM_ELEMENTOS;

    //cria a matriz
    for (i=0; i<tam; i++){
        for (j=0; j<tam; j++){
            M[i][j] = 1.0 / (i + j + 1.0);
		}
	}
    //cria o vetor fonte
    for (i=0; i<tam; i++){
        B[i] = 1.0 / (i + tam + 1.0);
	}

}

//Definicao da funcao solucao exata
double F(double x, double c1, double c2, double E){
    return c1 * exp(-1.0 * x / sqrt (E)) +  c2 * exp(x / sqrt (E)) + 1.0;
}
       
//##### codigo da questao 2 da lista 1 que gera a matriz do problema de valor de contorno #####

void gerarMatriz(double E, double M[NUM_ELEMENTOS][NUM_ELEMENTOS], double B[NUM_ELEMENTOS], double solucao_exata[NUM_ELEMENTOS], double tempo[NUM_ELEMENTOS]){

	int tam = NUM_ELEMENTOS;

    //Coeficientes
    double a = 0;
    double b = 1.0;
    int Nel = NUM_ELEMENTOS;
    double h = (b-a)/double(Nel);
    double c2 = (exp((-1.0)/sqrt(E)) - 1.0) / (exp(1.0/sqrt(E)) - exp((-1.0)/sqrt(E)));
    double c1 = - 1 - c2;
    
    double h2 = pow(h, 2);
    double cdp = 2.0 * E + h2;
    double dh = 0.0;
	int i, j;


	for (i=1; i<=tam; i++){
        dh = h * i;
        tempo[i-1]=dh;
        solucao_exata[i-1] = F(dh, c1, c2, E);
		//printf("i = %d, dh = %f\n", i, dh);
	}
/*
	for (i=0; i<num_particoes; i++){
		//printf("tempo[%d] = %f\n", i, tempo[i]);
		printf("solucao_exata[%d] = %f\n", i, solucao_exata[i]);
	}

	return;
	*/

	for (i=0; i<tam; i++){
		for (j=0; j<tam; j++){
            if(i==j){
                //preenche diagonal principal
                M[i][j] = cdp;
                B[i] = h2;
            } else if(i==j+1){
                //preenche diagonal inferior
                M[i][j] = E * -1.0;
            } else if(i==j-1){
                //preenche diagonal superior
                M[i][j] = E * -1.0;
            } else {
                //preenche o resto com zeros
                M[i][j] = 0;
			}
		}
	}
    
}

void retroSubstituicao(double M[NUM_ELEMENTOS][NUM_ELEMENTOS], double X[NUM_ELEMENTOS], double B[NUM_ELEMENTOS]){
    
	int n, i, j;
	int tam = NUM_ELEMENTOS;

    //checa se e superior
    bool isSuperior = (M[tam-1][0] == 0);

    //percorre as linhas da matriz (ta funcionando)
    for (n=0; n<tam; n++){
        
        //define o i se for superior ou inferior
        if(isSuperior){
            i = tam - n - 1;
		} else { 
            i = n;
		}
        
        //valor do vetor solucao da linha correspondente
        double temp = B[i];
        
        //soma os valores exceto o valor do pivo
        for (j=0; j<tam; j++){
            if(j != i){
                temp += X[j] * M[i][j] * -1;
			}
		}
       
        //interrompe se der divisao por zero
        double div =  M[i][i];
        if(div == 0){
            return;
		}
            
        //calcula a solucao da linha, usando a soma e o valor do pivo
        X[i] = temp / M[i][i];
	}
}

void gauss(double M[NUM_ELEMENTOS][NUM_ELEMENTOS], double X[NUM_ELEMENTOS], double B[NUM_ELEMENTOS]){
    
	int k, i, j, t;
	int tam = NUM_ELEMENTOS;
	double mult;

    //percorre os pivos da matriz zerando as linhas abaixo
    for (k=0; k<tam; k++){
        t = k + 1;
        
        //percorre os elementos abaixo do pivo para extrair o multiplicador
		for (i=t; i<tam; i++){
            mult = M[i][k] / M[k][k];
            
            //usa o multiplicador pra zerar o elemento da linha
			for (j=t; j<tam; j++){
                M[i][j] = M[i][j] - mult * M[k][j];
			}
            
            //usa o multiplicador pra mudar o elemento no vetor fonte
            B[i] = B[i] - mult * B[k];
		}
	}
    //agora que tem uma matriz diagonal superior, usa retroSubstituicao
    retroSubstituicao(M, X, B);
}

void gaussComPivoteamento(double M[NUM_ELEMENTOS][NUM_ELEMENTOS], double X[NUM_ELEMENTOS], double B[NUM_ELEMENTOS]){   
	int tam = NUM_ELEMENTOS;
	int a, b, i, j;
	
	//loop do pivoteamento
	for (a=0; a<tam; a++){

		//armazena pivo em modulo
		double maior_elemento = abs(M[a][a]);
		int linha_a = a;
		int linha_b = a;

		//percorre as colunas a procura do maior valor
		for (b=a; b<tam; b++){
			double elemento = abs(M[b][a]);
			if(elemento > maior_elemento){
				maior_elemento = elemento;
				linha_b = b;
			}
		}
		
		//se os indices das linhas a e b forem diferentes, faz a troca de linha
		if(linha_a != linha_b){
			trocarLinhas(M, B, linha_a, linha_b);
		}

		//percorre os elementos abaixo do pivo para extrair o multiplicador
		for (i=a+1; i<tam; i++){
			double mult = M[i][a] / M[a][a];
			
			//usa o multiplicador pra zerar o elemento da linha
			for (i=a; i<tam; i++){
				M[i][j] = M[i][j] - mult * M[a][j];
			}
			
			//usa o multiplicador pra mudar o elemento no vetor fonte
			B[i] = B[i] - mult * B[a];
		}
	}

	//agora que tem uma matriz diagonal superior, usa retroSubstituicao
	retroSubstituicao(M, X, B);
}

int main() {

	// cria variaveis
	double A[NUM_ELEMENTOS][NUM_ELEMENTOS];
	double B[NUM_ELEMENTOS];
	double X[NUM_ELEMENTOS];
	double XE[NUM_ELEMENTOS];
	double T[NUM_ELEMENTOS];

	//gerar a matriz
	gerarMatriz(0.01, A, B, XE, T);

	//imprimeMatriz(A);
	//return 0;

    //marca inicio da execucao

	gauss(A, X, B);

	//marca fim de execucao


	//gera o arquivo solucao pra salvar
	string arquivo =  "0|0|0\n";
	int i;
	for (i=0; i<NUM_ELEMENTOS; i++){
		arquivo = arquivo + to_string(T[i]) + "|" + to_string(XE[i]) + "|" + to_string(X[i]) + "\n";
	}
	ofstream out; // out é uma variavel.
	out.open("TesteData.txt");
	out<<arquivo;  // saida de uma variavel
	out.close(); // nã oesqueça de fechar...
	return 0;
}