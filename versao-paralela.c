#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define MESTRE 0

// DADOS COMPARTILHADOS
int m1[SIZE][SIZE], m2[SIZE][SIZE], mres[SIZE][SIZE];
int l1, c1, l2, c2, lres, cres;
int tamanhoMatriz = SIZE * SIZE;

int main(int argc, char *argv[])
{
  int id, p, hostsize;
  char hostname[MPI_MAX_PROCESSOR_NAME];
  double elapsed_time;
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  MPI_Get_processor_name(hostname, &hostsize);

  // INICIALIZA OS ARRAYS A SEREM MULTIPLICADOS.
  l1 = c1 = SIZE;
  l2 = c2 = SIZE;
  if (c1 != l2)
  {
    fprintf(stderr, "Impossivel multiplicar matrizes: parametros invalidos.\n");
    return 1;
  }
  lres = l1;
  cres = c2;

  // MESTRE.
  if (id == MESTRE)
  {
    // Popula as matrizes.
    int k = 1;
    for (int i = 0; i < SIZE; i++)
    {
      for (int j = 0; j < SIZE; j++)
      {
        if (k % 2 == 0)
          m1[i][j] = -k;
        else
          m1[i][j] = k;
      }
      k++;
    }
    k = 1;
    for (int j = 0; j < SIZE; j++)
    {
      for (int i = 0; i < SIZE; i++)
      {
        if (k % 2 == 0)
          m2[i][j] = -k;
        else
          m2[i][j] = k;
      }
      k++;
    }

    // PREPARA PARA MEDIR TEMPO
    // Comeca a contar depois que as matrizes estao populadas.
    elapsed_time = -MPI_Wtime();

    // BROADCAST PARA ENVIAR A MATRIZ 2 INTEIRA PARA TODOS OS ESCRAVOS
    MPI_Bcast(m2, tamanhoMatriz, MPI_INT, MESTRE, MPI_COMM_WORLD);

    // Loop para enviar partes da matriz 1 para os escravos.
    for (int i = 0; i < p - 1; ++i)
    {
      int idEscravo = i + 1;
      int idUltimoEscravo = p - 2;

      int qtdLinhas = SIZE / (p - 1);
      int numLinha = i * qtdLinhas;
      if (i == idUltimoEscravo)
      {
        // Se for o ultimo escravo, envia as linhas restantes.
        qtdLinhas += SIZE % (p - 1);
      }

      // SEND-NUMERO DA LINHA A COMECAR
      MPI_Send(&numLinha, 1, MPI_INT, idEscravo, 0, MPI_COMM_WORLD);

      // SEND-QUANTIDADE DE LINHAS A PROCESSAR
      MPI_Send(&qtdLinhas, 1, MPI_INT, idEscravo, 0, MPI_COMM_WORLD);

      // SEND-LINHAS M1
      MPI_Send(&m1[numLinha][0], qtdLinhas * SIZE, MPI_INT, idEscravo, 0, MPI_COMM_WORLD);
    }

    for (int i = 0; i < p - 1; ++i)
    {
      int numLinha, qtdLinhas;

      // RECEIVE-NUMERO DA LINHA A COMECAR
      MPI_Recv(&numLinha, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
      int source = status.MPI_SOURCE;

      // RECEIVE-QUANTIDADE DE LINHAS A PROCESSAR
      MPI_Recv(&qtdLinhas, 1, MPI_INT, source, 0, MPI_COMM_WORLD, &status);
      
      // RECEIVE-LINHAS M RESULTANTE
      MPI_Recv(&mres[numLinha][0], qtdLinhas * SIZE, MPI_INT, source, 0, MPI_COMM_WORLD, &status);
    }

    // OBTEM O TEMPO
    // Para de contar depois que a matriz resultante esta pronta.
    elapsed_time += MPI_Wtime();
    printf("%d - %lf\n", SIZE, elapsed_time);

    // VERIFICA SE O RESULTADO DA MULTIPLICACAO ESTA CORRETO
    for (int i = 0; i < SIZE; i++)
    {
      int k = SIZE * (i + 1);
      for (int j = 0; j < SIZE; j++)
      {
        int k_col = k * (j + 1);
        if (i % 2 == 0)
        {
          if (j % 2 == 0)
          {
            if (mres[i][j] != k_col)
              printf("\nERRO validando matriz resultante: %d != %d", mres[i][j], k_col);
              return 1;
          }
          else
          {
            if (mres[i][j] != -k_col)
              printf("\nERRO validando matriz resultante: %d != %d", mres[i][j], -k_col);
              return 1;
          }
        }
        else
        {
          if (j % 2 == 0)
          {
            if (mres[i][j] != -k_col)
              printf("\nERRO validando matriz resultante: %d != %d", mres[i][j], -k_col);
              return 1;
          }
          else
          {
            if (mres[i][j] != k_col)
              printf("\nERRO validando matriz resultante: %d != %d", mres[i][j], k_col);
              return 1;
          }
        }
      }
    }
  }
  // ESCRAVO
  else
  {
    // Recebe a matriz 2 inteira.
    MPI_Bcast(m2, tamanhoMatriz, MPI_INT, MESTRE, MPI_COMM_WORLD);

    int qtdLinhas, numLinha;

    // RECEIVE-NUMERO DA LINHA A COMECAR
    MPI_Recv(&numLinha, 1, MPI_INT, MESTRE, 0, MPI_COMM_WORLD, &status);

    // RECEIVE-QUANTIDADE DE LINHAS A PROCESSAR
    MPI_Recv(&qtdLinhas, 1, MPI_INT, MESTRE, 0, MPI_COMM_WORLD, &status);

    // RECEIVE-LINHAS M1
    MPI_Recv(&m1[numLinha][0], qtdLinhas * SIZE, MPI_INT, MESTRE, 0, MPI_COMM_WORLD, &status);

    // REALIZA A MULTIPLICACAO
    // Utiliza OpenMP para paralelizar o loop.
#pragma omp parallel for
    for (int i = numLinha; i < numLinha + qtdLinhas; i++)
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

    // SEND-NUMERO DA LINHA A COMECAR
    MPI_Send(&numLinha, 1, MPI_INT, MESTRE, 0, MPI_COMM_WORLD);

    // SEND-QUANTIDADE DE LINHAS A PROCESSAR
    MPI_Send(&qtdLinhas, 1, MPI_INT, MESTRE, 0, MPI_COMM_WORLD);

    // SEND-LINHAS M RESULTANTE
    MPI_Send(&mres[numLinha][0], qtdLinhas * SIZE, MPI_INT, MESTRE, 0, MPI_COMM_WORLD);
  }

  MPI_Finalize();
  return 0;
}