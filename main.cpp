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
  int bquan, qquan = 0;
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
  bquan *= m;
  quan *= (n * m);
  int re = k%p;
  if (my_rank < re)
  {
      quan += (n * m);
      bquan += m;
      qquan++;
  }
  if(my_rank == re && l)
  {
    bquan += l;
    qquan++;
    quan+= (l * n);
  }
  double *a = new double[quan + 3 * bquan + n];
  if (!a)
    {
      if (my_rank == 0)
        printf ("Not enough memory\n");
      return 0;
    }
  double *b = a + quan, *r = b + bquan, *x = r + bquan, *x1 = x + bquan;
  if (namea)
    {
      res = read_matrix (a, n, m, namea, my_rank, p);
      if (res < 0)
        {
          delete []a;
          if (my_rank == 0)
            {
              if (res <= -1 && res >= -1)
                {
                  printf ("Couldn't open file %s\n", namea);
                }
              else
                {
                  printf ("Couldn't read from file %s\n", namea);
                }
            }
          MPI_Finalize ();
          return 0;
        }
      else
        {
          norm = res;
        } 
    }
  else
    norm = init_matrix (a, n, m, my_rank, p);
  print_matrix (a, n, m, my_rank, p);
  if (my_rank == 0)
    printf ("%.2f\n", norm);
  build_b (x, n, a, m, my_rank, p);
  build_x (b, n, m, my_rank, p);
  build_r (r, n, m, my_rank, p);
  MPI_Barrier (MPI_COMM_WORLD);
  /*for (int i = 0; i < p; i++)
      {
          if (i % p == my_rank)
          {
              for (int j = 0; j < bquan; j++)
                  printf ("%.2f\n", x[j]);
          }
      }*/
  res = solve (a, x, b, n, m, my_rank, p, norm);
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
  if (namea)
    {
      res = read_matrix (a, n, m, namea, my_rank, p);
      if (res < 0)
        {
          delete []a;
          if (my_rank == 0)
            {
              if (res <= -1 && res >= -1)
                {
                  printf ("Couldn't open file %s\n", namea);
                }
              else
                {
                  printf ("Couldn't read from file %s\n", namea);
                }
            }
          MPI_Finalize ();
          return 0;
        }
      else
        {
          norm = res;
        }
    }
  else
    norm = init_matrix (a, n, m, my_rank, p);
  /*for (int i = 0; i < p; i++)
    {
        if (i % p == my_rank)
        {
            for (int j = 0; j < bquan; j++)
                printf ("%.2f\n", b[j]);
        }
    }*/
  build_b (b, n, a, m, my_rank, p);
  for (int i = 0; i < p; i++)
  {
      if (i % p == my_rank)
      {
          for (int j = 0; j < bquan; j++)
              printf ("%.2f\n", x[j]);
      }
  }
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
