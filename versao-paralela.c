#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

// DADOS COMPARTILHADOS
int m1[SIZE][SIZE],m2[SIZE][SIZE],mres[SIZE][SIZE];
int l1, c1, l2, c2, lres, cres;

int main(int argc, char *argv[]) {
  int i, j, k, id, p;
  double elapsed_time;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  MPI_Status status;
  /* Só faz sentido na versão sequencial
  if (id != 0) {
      MPI_Finalize();
      exit(0);
  }
  */

  // INICIALIZA OS ARRAYS A SEREM MULTIPLICADOS.
  l1 = c1 = SIZE;
  l2 = c2 = SIZE;
  if (c1 != l2) {
    fprintf(stderr, "Impossivel multiplicar matrizes: parametros invalidos.\n");
    return 1;
  }
  lres = l1;
  cres = c2;

  //Verificar se é 0 para saber se é o mestre
  if (id == 0) {
    //Se for o mestre, precisa fazer o trabalho inicial de preencher as matrizes
    int k = 1;
    for (i=0; i<SIZE; i++) {
      for (j=0 ; j<SIZE; j++) {
        if (k%2==0)
          m1[i][j] = -k;
        else
          m1[i][j] = k;
      }
      k++;
    }
    k = 1;
    for (j=0 ; j<SIZE; j++) {
        for (i=0 ; i<SIZE; i++) {
	        if (k%2==0)
               m2[i][j] = -k;
	        else
               m2[i][j] = k;
        }
        k++;
    }

    // PREPARA PARA MEDIR TEMPO
    elapsed_time = -MPI_Wtime();

    //Todos os escravos precisam ter uma das matrizes completas, escolhemos mandar a M2
    MPI_Bcast(m2, pow(SIZE,2), MPI_INT, 0, MPI_COMM_WORLD);

    //Enviar parte da m1 para os escravos
    for (i = 1; i < p; ++i) {
      int idEscravo = i; 

      int numeroDeLinhasPorProcesso = SIZE / (p - 1);
      int linhaInicialDoProcesso = (i-1) * numeroDeLinhasPorProcesso;
      if (i == (p - 1)) { //Verificar se é o último processo a ser enviado os dados
        //Mandar somente o "resto", linhas restantes para o último
        numeroDeLinhasPorProcesso += SIZE % (p - 1);
      }

      // Enviar a linha que o processo deve começar a calcular, a quantidade de linhas e as linhas da matriz
      MPI_Send(&linhaInicialDoProcesso, 1, MPI_INT, idEscravo, 0, MPI_COMM_WORLD);
      MPI_Send(&numeroDeLinhasPorProcesso, 1, MPI_INT, idEscravo, 0, MPI_COMM_WORLD);
      MPI_Send(&m1[linhaInicialDoProcesso][0], numeroDeLinhasPorProcesso * SIZE, MPI_INT, idEscravo, 0, MPI_COMM_WORLD);
    }

    //For para pegar os retornos após o cálculos em cada escravo
    for (int i = 1; i < p; ++i) {
      int numeroDeLinhasPorProcesso;
      int linhaInicialDoProcesso;

      //Usando o MPI_Recv para receber primeiro o numero da linha do intervalo que foi calculado
      //Segundo a quantidade de linhas que vai ser fornecida no retorno e terceiro os valores da matriz resultante
      MPI_Recv(&linhaInicialDoProcesso, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(&numeroDeLinhasPorProcesso, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(&mres[linhaInicialDoProcesso][0], numeroDeLinhasPorProcesso * SIZE, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &status);
    }

    // OBTEM O TEMPO
    elapsed_time += MPI_Wtime();
    printf("%d - %lf\n", SIZE, elapsed_time);

    // VERIFICA SE O RESULTADO DA MULTIPLICACAO ESTA CORRETO
    for (i=0 ; i<SIZE; i++) {
      k = SIZE*(i+1);
      for (j=0 ; j<SIZE; j++) {
          int k_col = k*(j+1);
          if (i % 2 ==0) {
             if (j % 2 == 0) {
                if (mres[i][j]!=k_col)
                   return 1;
             }
             else {
                if (mres[i][j]!=-k_col)
                   return 1;
             }
          }
          else {
             if (j % 2 == 0) {
                if (mres[i][j]!=-k_col)
                   return 1;
             }
             else {
                if (mres[i][j]!=k_col)
                   return 1;
             }
          }
      } 
    }
  }
  // Se for diferente de 0, significa que o processo é um escravo
    else
  {
    // Recebe a matriz 2 inteira.
    MPI_Bcast(m2, (SIZE * SIZE), MPI_INT, 0, MPI_COMM_WORLD);

    int numeroDeLinhasPorProcesso;
    int linhaInicialDoProcesso;

    MPI_Recv(&linhaInicialDoProcesso, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    MPI_Recv(&numeroDeLinhasPorProcesso, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    MPI_Recv(&m1[linhaInicialDoProcesso][0], numeroDeLinhasPorProcesso * SIZE, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

    // REALIZA A MULTIPLICACAO
    for (int i = linhaInicialDoProcesso; i < linhaInicialDoProcesso + numeroDeLinhasPorProcesso; i++)
    {
      for (int j = 0; j < cres; j++)
      {
        mres[i][j] = 0;
        for (int k = 0; k < c1; k++)
        {
          mres[i][j] += m1[i][k] * m2[k][j];
        }
      }
    }

    MPI_Send(&numLinha, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&numeroDeLinhasPorProcesso, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&mres[linhaInicialDoProcesso][0], numeroDeLinhasPorProcesso * SIZE, MPI_INT, 0, 0, MPI_COMM_WORLD);
  }
  MPI_Finalize();
  return 0;
}