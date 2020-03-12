#include "matrixes.h"
#include "time.h"
#include <sys/sysinfo.h>

int
main (int argc, char *argv[])
{
  int n, m;
  double res = 0;
  char *namea = 0;
  double norm = 0;
  int my_rank; //номер текущего процесса
  int p; // количество процессов
  int k,
      l; // n / m
  int quan; // колво требуемой памяти для текущего процесса(на матрицу)
  int bquan, qquan = 0, block_quan = 0;
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
  block_quan = 3 * m * m + 2 + n * m;
  double *a = new double[quan + 2 * bquan + n + block_quan];
  if (!a)
    {
      if (my_rank == 0)
        printf ("Not enough memory\n");
      MPI_Finalize ();
      return 0;
    }
  double *x = a + quan, *b = x + bquan, *x1 = b + bquan, *block = x1 + n;
  norm = init_forms (namea, a, n, m, my_rank, p);
  if (norm <= -1 && norm >= -1)
    {
      delete []a;
      MPI_Finalize ();
      return 0;
    }
  print_matrix (a, n, m, my_rank, p);
  if (my_rank == 0)
    printf ("%.2f\n", norm);
  build_b (x, n, a, m, my_rank, p);
  MPI_Barrier (MPI_COMM_WORLD);
  double t = MPI_Wtime ();
  res = solve (a, x, n, m, my_rank, p, norm, block);
  MPI_Barrier (MPI_COMM_WORLD);
  t = MPI_Wtime () - t;
  if (res < 0)
    {
      if (my_rank == 0)
        printf ("Solution cannot be found!\n");
      delete []a;
      MPI_Finalize ();
      return 0;
    }
  MPI_Barrier (MPI_COMM_WORLD);
  init_forms (namea, a, n, m, my_rank, p);
  build_b (b, n, a, m, my_rank, p);
  MPI_Barrier (MPI_COMM_WORLD);
  double Residual = AX_B_MPI (a, x1, n, my_rank, p, m, b, x);
  MPI_Barrier (MPI_COMM_WORLD);
  if (my_rank == 0)
    printf ("Residual = %e Elapsed = %.2f n = %d m = %d p = %d\n", Residual, t, n, m, p);
  //print_matrix (a, n, m, my_rank, p);
  fflush (stdout);
  fflush (stderr);
  delete []a;
  MPI_Finalize ();
  return 0;
}
