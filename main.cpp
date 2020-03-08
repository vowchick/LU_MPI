#include "matrixes.h"
#include "time.h"

int
main (int argc, char *argv[])
{
  int n, m;
  double res;
  char *namea = 0;
  double norm = 0;
  int my_rank; //номер текущего процесса
  int p; // количество процессов
  int k,
      l; // n / m
  int quan; // колво требуемой памяти для текущего процесса(на матрицу)
  int bquan, qquan = 0, block_quan = 0, tot_quan= 0;
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &p);
  if (argc < 3 || argc > 4 || ((n = atoi (argv[1])) <= 0) || ((m = atoi (argv[2])) <= 0))
    {
      if (my_rank == 0)
        printf ("usage: %s n (>0) m(>0) [namea] \n", argv[0]);
      MPI_Finalize ();
      return 1;
    }
  if (n < m)
    {
      if (my_rank == 0)
        printf ("n should be such as: n >= m\n");
      MPI_Finalize ();
      return 0;
    }
  if (m % 3 != 0)
    {
      if (my_rank == 0)
        printf ("Pleas input m so m (mod 3) = 0\n");
      MPI_Finalize ();
      return 0;
    }
  if (argc == 4)
    namea = argv[3];
  k = n / m;
  l = n % m;
  quan = k / p;
  qquan = k / p;
  bquan = k / p;
  tot_quan = k / p;
  tot_quan *= (m * m);
  bquan *= m;
  quan *= (n * m);
  int re = k%p;
  if (my_rank < re)
  {
      quan += (n * m);
      bquan += m;
      qquan++;
      tot_quan += (m * m);
  }
  if(my_rank == re && l)
  {
    bquan += l;
    qquan++;
    quan+= (l * n);
    tot_quan += (m * l);
  }
  block_quan = 2 * m * m + 1 + tot_quan + n * m;
  double *a = new double[quan + 2 * bquan + n + block_quan];
  if (!a)
    {
      if (my_rank == 0)
        printf ("Not enough memory\n");
      return 0;
    }
  double *x = a + quan, *b = x + bquan, *x1 = b + bquan, *block = x1 + n;
  norm = init_forms (namea, a, n, m, my_rank, p);
  if (norm <= -1 && norm >= -1)
    {
      MPI_Finalize ();
      return 0;
    }
  print_matrix (a, n, m, my_rank, p);
  if (my_rank == 0)
    printf ("%.2f\n", norm);
  build_b (x, n, a, m, my_rank, p);
  MPI_Barrier (MPI_COMM_WORLD);
  /*for (int i = 0; i < p; i++)
      {
          if (i % p == my_rank)
          {
              for (int j = 0; j < bquan; j++)
                  printf ("%.2f\n", x[j]);
          }
      }*/
  res = solve (a, x, n, m, my_rank, p, norm, block);
   MPI_Barrier (MPI_COMM_WORLD);
  if (res < 0)
    {
      if (my_rank == 0)
        printf ("Some awfull mistake occured!\n");
      delete []a;
      MPI_Finalize ();
      return 0;
    }
  MPI_Barrier (MPI_COMM_WORLD);
  init_forms (namea, a, n, m, my_rank, p);
  build_b (b, n, a, m, my_rank, p);
  /*for (int i = 0; i < p; i++)
  {
      if (i % p == my_rank)
      {
          for (int j = 0; j < bquan; j++)
              printf ("%.2f\n", x[j]);
      }
  }*/
  MPI_Barrier (MPI_COMM_WORLD);
  double Residual = AX_B_MPI (a, x1, n, my_rank, p, m, b, x);
  MPI_Barrier (MPI_COMM_WORLD);
  if (my_rank == 0)
    printf ("Residual = %e\n", Residual);
  print_matrix (a, n, m, my_rank, p);
  fflush (stdout);
  fflush (stderr);
  delete []a;
  MPI_Finalize ();
  return 0;
}
