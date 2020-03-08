#include "matrixes.h"
double
read_matrix (double *a, int n, int m, const char *name, int my_rank, int p)
{
  int k = n / m, l = n - k * m;
  int i, j, s, q, loc_q = 0;
  double *buf;
  double sum = 0;
  double max[1];
  double rmax[1];
  max[0] = 0;
  int err[1];
  int rerr[1];
  int flag = 0;
  err[0] = 0;
  FILE *fp = 0;
  MPI_Status status;
  int tag = 0;
  buf = new double [n * m];
  if (my_rank == 0)
    {
      fp = fopen (name, "r");
      if (!fp)
        {
          err[0] = -1;
          flag = -1;
        }
    }
  MPI_Allreduce(err, rerr, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  for (q = 0; q < k && rerr[0] == 0; q++)
    {
      if (my_rank == 0)
        {
          for (s = 0; s < m && err[0] == 0; s++)
            {
              for (j = 0; j < k && err[0] == 0; j++)
                {
                  for (i = 0; i < m && err[0] == 0; i++)
                    {
                      if (q % p == 0)
                        {
                          if (fscanf (fp, "%lf", a + i + j * m * m + s * m + loc_q * m * n) != 1)
                            {
                              err[0] = -2;
                              break;
                            }
                          sum += fabs (a[i + j * m * m + s * m + loc_q * m * n]);
                        }
                      else
                        {
                          if (fscanf (fp, "%lf", buf + i + j * m * m + s * m) != 1)
                            {
                              err[0] = -2;
                              break;
                            }
                          sum += fabs (buf[i + j * m * m + s * m]);
                        }
                    }
                }
              for (i = 0; i < l && err[0] == 0; i++)
                    {
                      if (q % p == 0)
                        {
                          if (fscanf (fp, "%lf", a + i + k * m * m + s * l + loc_q * m * n) != 1)
                            {
                              err[0] = -2;
                              break;
                            }
                        }
                      else
                        {
                          if (fscanf (fp, "%lf", buf + i + k * m * m + s * l) != 1)
                            {
                              err[0] = -2;
                              break;
                            }
                        }
                    }
              if (max[0] < sum)
                max[0] = sum;
              sum = 0;
            }
          if (q % p != 0)
            MPI_Send (buf, n * m, MPI_DOUBLE,
                      q % p, tag, MPI_COMM_WORLD);
          else
            loc_q++;
        }
        else
          {
            if (my_rank == q % p)
              {
                MPI_Recv (a + loc_q * n * m, n * m, MPI_DOUBLE,
                          0, tag, MPI_COMM_WORLD, &status);
                loc_q++;
              }
          }
        MPI_Allreduce(err, rerr, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    }
  MPI_Barrier (MPI_COMM_WORLD);
    if (rerr[0] == 0)
    {
      if (my_rank == 0)
          {
            for (s = 0; s < l && err[0] == 0; s++)
              {
                for (j = 0; j < k && err[0] == 0; j++)
                  {
                    for (i = 0; i < m && err[0] == 0; i++)
                      {
                        if (k % p == 0)
                        {
                          if (fscanf (fp, "%lf", a + i + j * m * l + s * m + loc_q * m * n) != 1)
                            {
                              err[0] = -2;
                              break;
                            }
                          sum += fabs (a[i + j * m * l + s * m + loc_q * m * n]);
                        }
                        else
                          {
                            if (fscanf (fp, "%lf", buf + i + j * m * l + s * m) != 1)
                              {
                                err[0] = -2;
                                break;
                              }
                            sum += fabs (buf[i + j * m * l + s * m]);
                          }
                      }
                  }
                for (i = 0; i < l && rerr[0] == 0; i++)
                  {
                    if (k % p == 0)
                      {
                        if (fscanf (fp, "%lf", a + i + k * m * l + s * l + loc_q * m * n) != 1)
                          {
                            err[0] = -2;
                            break;
                          }
                        sum += fabs (a[i + k * m * l + s * l + loc_q * m * n]);
                      }
                    else
                      {
                        if (fscanf (fp, "%lf", buf + i + k * m * l + s * l) != 1)
                          {
                            err[0] = -2;
                            break;
                          }
                        sum += fabs (buf[i + k * m * l + s * l]);
                      }
                  }
                if (max[0] < sum)
                  max[0] = sum;
                sum = 0;
              }
            if (k % p != 0)
              MPI_Send (buf, n * l, MPI_DOUBLE,
                        k % p, tag, MPI_COMM_WORLD);
            else
              loc_q++;
          }
          else
            {
              if (my_rank == k % p)
                {
                  MPI_Recv (a + loc_q * n * m, n * l, MPI_DOUBLE,
                            0, tag, MPI_COMM_WORLD, &status);
                  loc_q++;
                }
            }
    }
  if (my_rank == 0 && flag == 0)
    fclose (fp);
  delete []buf;
  MPI_Allreduce(err, rerr, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce (max, rmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  if (rerr[0] < 0)
    return rerr[0];
  return rmax[0];
}
double
init_forms (char *namea, double *a, int n, int m, int my_rank, int p)
{
    double norm = 0, res = 0;
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
            return -1;
          }
        else
          {
            norm = res;
          }
      }
    else
      norm = init_matrix (a, n, m, my_rank, p);
  return norm;
}
void
print_matrix (double *a, int n, int m, int my_rank, int p)
{
  int q, k = n / m, loc_q = 0, l = n % m, quan;
  int tag = 0;
  double *buf;
  if (my_rank == 0)
    buf = new double[m * n];
  MPI_Status status;
  for (q = 0; q < k && q * m < max_print; q++)
    {
      if (my_rank == q % p && q % p != 0)
        {
          MPI_Send (a + loc_q * n * m, n * m, MPI_DOUBLE,
                    0, tag, MPI_COMM_WORLD);
          loc_q++;
        }
      else if (my_rank == 0)
        {
          if (q % p != 0)
            {
              MPI_Recv (buf, n * m, MPI_DOUBLE, q % p, tag,
                        MPI_COMM_WORLD, &status);
              if (q * m + m > max_print)
                {
                  quan = max_print - q * m;
                }
              else
                quan = m;
              print_block_string (buf, n, m, m, quan);
            }
          else
            {
              if (q * m + m > max_print)
                {
                  quan = max_print - q * m;
                }
              else
                quan = m;
              print_block_string (a + loc_q * n * m, n, m, m, quan);
              loc_q++;
            }
        }
    }
    if (l != 0 && k % p == my_rank && my_rank != 0)
      {
        MPI_Send (a + loc_q * n * m, n * l, MPI_DOUBLE,
                    0, tag, MPI_COMM_WORLD);
          loc_q++;
      }
    else if (l != 0 && my_rank == 0)
      {
        if (k % p != 0)
            {
              MPI_Recv (buf, n * l, MPI_DOUBLE, k % p, tag,
                        MPI_COMM_WORLD, &status);
              if (k * m + m > max_print)
                {
                  quan = max_print - k * m;
                }
              else
                quan = l;
              print_block_string (buf, n, m, l, quan);
            }
          else
            {
              if (k * m + m > max_print)
                {
                  quan = max_print - k * m;
                }
              else
                quan = l;
              print_block_string (a + loc_q * n * m, n, m, l, quan);
              loc_q++;
            }
      }
  if (my_rank == 0)
    delete []buf;
}
void
print_block_string (double *a, int n, int m, int m_l, int quan)
{
   int q, i, j, k = n / m, l = n % m;
   for (q = 0; q < quan; q++)
        {
           for (i = 0; i < k; i++)
            {
                for (j = 0; j < m && j + i * m < max_print; j++)
                    printf ("%.2f ", a[j + i * m * m_l + q * m]);

            }
           for (j = 0; j < l && j + k * m < max_print; j++)
             {
                printf ("%.2f ", a[j + k * m * m_l + q * l]);
             }
            printf ("\n");

        }
}
double
init_matrix (double *a, int n, int m, int my_rank, int p)
{
  int i, j, q, s, k = n / m, l = n - k * m, loc_q = 0;
  int tag = 0;
  double sum = 0;
  double max[1];
  double rmax[1];
  max[0] = 0;
  MPI_Status status;
  double *buf = new double [n * m];
  for (q = 0; q < k; q++)
    {
      if (my_rank == 0)
        {
          for (s = 0; s < m; s++)
            {
              for (j = 0; j < k; j++)
                {
                  for (i = 0; i < m; i++)
                    {
                      if (q % p == 0)
                        {
                          a[i + j * m * m + s * m + loc_q * m * n] = f1 (q * m + s, i + j * m, n);
                          sum += fabs(a[i + j * m * m + s * m + loc_q * m * n]);
                        }
                      else
                        {
                          buf[i + j * m * m + s * m] = f1 (q * m + s, i + j * m, n);
                          sum += fabs(buf[i + j * m * m + s * m]);
                        }
                    }
                }
              for (i = 0; i < l; i++)
                    {
                      if (q % p == 0)
                        {
                          a[i + k * m * m + s * l + loc_q * m * n] = f1 (q * m + s, k * m + i, n);
                          sum += fabs (a[i + k * m * m + s * l + loc_q * m * n]);
                        }
                      else
                        {
                          buf[i + k * m * m + s * l] = f1 (q * m + s, k * m + i, n);
                          sum += fabs(buf[i + k * m * m + s * l]);
                        }
                    }
              if (max[0] < sum)
                max[0] = sum;
              sum = 0;
            }
          if (q % p != 0)
            MPI_Send (buf, n * m, MPI_DOUBLE,
                      q % p, tag, MPI_COMM_WORLD);
          else
            loc_q++;
        }
        else
          {
            if (my_rank == q % p)
              {
                MPI_Recv (a + loc_q * n * m, n * m, MPI_DOUBLE,
                          0, tag, MPI_COMM_WORLD, &status);
                loc_q++;
              }
          }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (my_rank == 0)
        {
          for (s = 0; s < l; s++)
            {
              for (j = 0; j < k; j++)
                {
                  for (i = 0; i < m; i++)
                    {
                      if (k % p == 0)
                      {
                        a[i + j * m * l + s * m + loc_q * m * n] = f1 (k * m + s, i + j * m, n);
                        sum += fabs (a[i + j * m * l + s * m + loc_q * m * n]);
                      }
                      else
                        {
                          buf[i + j * m * l + s * m] = f1 (k * m + s, i + j * m, n);
                          sum += fabs (buf[i + j * m * l + s * m]);
                        }
                    }
                }
              for (i = 0; i < l; i++)
                {
                  if (k % p == 0)
                    {
                      a[i + k * m * l + s * l + loc_q * m * n] = f1 (k * m + s, k * m + i, n);
                      sum += fabs (a[i + k * m * l + s * l + loc_q * m * n]);
                    }
                  else
                    {
                      buf[i + k * m * l + s * l] = f1 (k * m + s, k * m + i, n);
                      sum += fabs (buf[i + k * m * l + s * l]);
                    }
                }
              if (max[0] < sum)
                max[0] = sum;
              sum = 0;
            }
          if (k % p != 0)
            MPI_Send (buf, n * l, MPI_DOUBLE,
                      k % p, tag, MPI_COMM_WORLD);
          else
            loc_q++;
        }
        else
          {
            if (my_rank == k % p)
              {
                MPI_Recv (a + loc_q * n * m, n * l, MPI_DOUBLE,
                          0, tag, MPI_COMM_WORLD, &status);
                loc_q++;
              }
          }
  delete []buf;
  MPI_Allreduce (max, rmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return rmax[0];
}
void
build_b (double *b, int n, double *a, int m, int my_rank, int p)
{
  int i, j, q, s, k = n / m, l = n - k * m, loc_q = 0;
  double sum = 0;
  for (q = 0; q < k; q++)
    {
      if (my_rank == q % p)
      {
        for (s = 0; s < m; s++)
          {
            for (j = 0; j < k; j++)
              {
                for (i = 0; i < m; i++)
                  {
                    sum += a[i + j * m * m + s * m + loc_q * m * n] * ((i + j * m) % 2);
                  }
              }
            for (i = 0; i < l; i++)
                  {
                    sum += a[i + k * m * m + s * l + loc_q * m * n] * ((i + k * m) % 2);
                  }
            b[loc_q * m + s] = sum;
            //printf ("%.2f\n", sum);
            sum = 0;
          }
        loc_q++;
      }
    }
    if (my_rank == k % p)
      {
        for (s = 0; s < l; s++)
          {
            for (j = 0; j < k; j++)
              {
                for (i = 0; i < m; i++)
                  {
                    sum += a[i + j * m * l + s * m + loc_q * m * n] * ((i + j * m) % 2);
                  }
              }
            for (i = 0; i < l; i++)
              {
                sum += a[i + k * m * l + s * l + loc_q * m * n] * ((i + k * m) % 2);
              }
            b[loc_q * m + s] = sum;
            sum = 0;
             //printf ("%.2f\n", sum);
          }
        loc_q++;
      }
}
double
f1 (int i, int j, int n)
{
  //return 1./(i + j + 1);
  if (i > j)
    return n - i;
  return n - j;
}
void
print_mistakes (char *name, int res)
{
  switch (res)
    {
      case -1:
        {
          printf ("Couldn't open file %s\n", name);
          break;
        }
      case -2:
        {
          printf ("Couldn't read from file %s\n", name);
          break;
        }
      default:
        printf ("Something else is wrong\n");
    }
}
void
build_r (double *x, int n, int m, int my_rank, int p)
{
  int k = n / m, l = n % m, loc_q = 0;
  for (int i = 0; i < k; i++)
    {
      if (my_rank == i % p)
        {
          for (int j = 0; j < m; j++)
            {
              x[loc_q * m + j] = ((i * m + j) % 2 == 0? 1:0);
            }
          loc_q++;
        }
    }
    if (my_rank == k % p)
      {
        for (int j = 0; j < l; j++)
          {
            x[loc_q * m + j] = ((k * m + j) % 2 == 0? 1:0);
          }
        loc_q++;
      }
}
void
build_x (double *x, int n, int m, int my_rank, int p)
{
  int k = n / m, l = n % m, loc_q = 0;
  for (int i = 0; i < k; i++)
    {
      if (my_rank == i % p)
        {
          for (int j = 0; j < m; j++)
            {
              x[loc_q * m + j] = 0.;
            }
          loc_q++;
        }
    }
    if (my_rank == k % p)
      {
        for (int j = 0; j < l; j++)
          {
            x[loc_q * m + j] = 0.;
          }
        loc_q++;
      }
}
