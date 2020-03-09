﻿#include "matrixes.h"
//#include "time.h"
//in case of a mistake should return -10
int
solve (double *a, double *x, int n, int m, int my_rank, int p, double norrm, double *block)
{
    int i, j, q, k = n / m, l = n % m;
    int str = 0, str2 = 0, str3 = 0;
    int quan = k / p, tot_quan = k / p, an_quan = k / p;
    tot_quan *= (m * m);
    an_quan *= m;
    int re = k%p;
    if (my_rank < re)
    {
        quan++;
        tot_quan += (m * m);
        an_quan += m;
    }
    if(my_rank == re && l)
    {
        quan++;
        an_quan += l;
        tot_quan += (m * l);
    }
    double *inv, *big_row;
    inv = block + m * m;
    big_row = inv + m * m + 1;
    for (i = 1; i <= k; i++)
      {
        for (q = 0; q < i; q++)
          {
            if (my_rank == q % p)
              {
                find_inv (q + 1, my_rank, p, n, m, l, quan, re, a, block, inv, norrm, m);
                put_blocks_into_string (a, big_row, block, q, m, n, quan, re, my_rank, p);
                put_inv(big_row, n * m, inv, m);
                MPI_Bcast (big_row, (n + m) * m + 1, MPI_DOUBLE, q % p, MPI_COMM_WORLD);
                if (big_row[(n + m) * m] < 0)
                  {
                    printf ("something\n");
                    return big_row[(n + m) * m];
                  }
              }
            else
              {
                  MPI_Bcast (big_row, (n + m) * m + 1, MPI_DOUBLE, q % p, MPI_COMM_WORLD);
                  if (big_row[(n + m) * m] < 0)
                    {
                      printf ("something\n");
                      return big_row[(n + m) * m];
                    }
              }
            if (i % p == my_rank)
              {
                if (i != k)
                {
                    str = get_bounds (i, q, n, m, l, quan, re, my_rank, p);
                    get_block(a, block, str, m, m);
                    str2 = n * m;
                    sq_prod (block, big_row + str2, a + str, m);
                    for (j = q + 1; j < k; j++)
                      {
                        str = get_bounds (i, q, n, m, l, quan, re, my_rank, p);
                        str2 = get_bounds_row (q, j, n, m);
                        sq_prod (a + str, big_row + str2, block, m);
                        str3 = get_bounds (i, j, n, m, l, quan, re, my_rank, p);
                        sum (a + str3, block, m);
                      }
                    if (l)
                      {
                        str = get_bounds (i, q, n, m, l, quan, re, my_rank, p);
                        str2 = get_bounds_row (q, k, n, m);
                        mult (a + str, m, m, l, big_row + str2, block);
                        str3 = get_bounds (i, k, n, m, l, quan, re, my_rank, p);
                        sum2 (a + str3, block, -1, m, l);
                      }
                }
                else
                  {
                    str = get_bounds (k, q, n, m, l, quan, re, my_rank, p);
                    get_block(a, block, str, l, m);
                    str2 = n * m;
                    mult (block, l, m, m, big_row + str2, a + str);
                    //sq_prod (block, big_row + str2, a + str, m);
                    for (j = q + 1; j < k; j++)
                      {
                        str = get_bounds (k, q, n, m, l, quan, re, my_rank, p);
                        str2 = get_bounds_row (q, j, n, m);
                        mult (a + str, l, m, m, big_row + str2, block);
                        //sq_prod (a + str, big_row + str2, block, m);
                        str3 = get_bounds (k, j, n, m, l, quan, re, my_rank, p);
                        sum2 (a + str3, block, -1, l, m);
                        //sum (a + str3, block, m);
                      }
                    if (l)
                      {
                        str = get_bounds (k, q, n, m, l, quan, re, my_rank, p);
                        str2 = get_bounds_row (q, k, n, m);
                        mult (a + str, l, m, l, big_row + str2, block);
                        str3 = get_bounds (k, k, n, m, l, quan, re, my_rank, p);
                        sum (a + str3, block, l);
                        //sum2 (a + str3, block, -1, m, l);
                      }
                  }
              }
          }
      }
    print_matrix(a, n, m, my_rank, p);
    MPI_Barrier (MPI_COMM_WORLD);
    Reverse (a, x, n, m, my_rank, p, norrm, re, quan, block);
    return 0;
}
int
Reverse (double *a, double *b,
         int n, int m, int my_rank, int p, double norrm, int re, int quan, double *block)
{
    int k = n / m, l = n % m, j = 0, str = 0, i = 0, str2 = 0;
    double *inv = block + m * m, *column = inv + m * m + 1, *buf = column + m;
    for (j = 0; j < k; j++)
    {
        if (my_rank == j % p)
        {
            str = get_bounds_vector (j, m, p);
            get_vect(b, buf, str, str + m);
            MPI_Bcast (buf, m, MPI_DOUBLE, my_rank, MPI_COMM_WORLD);
        }
        else
        {
            MPI_Bcast (buf, m, MPI_DOUBLE, j % p, MPI_COMM_WORLD);
        }
        for (i = j + 1; i < k; i++)
        {
            if (my_rank == i % p)
            {
                str = get_bounds_vector (i, m, p);
                str2 = get_bounds (i, j, n, m, l, quan, re, my_rank, p);
                mult (a + str2, m, m, 1, buf, column);
                v_sum (b + str, column, -1, m);
            }
        }
        if (l)
        {
            if (my_rank == k % p)
            {
                str2 = get_bounds (k, j, n, m, l, quan, re, my_rank, p);
                mult (a + str2, l, m, 1, buf, column);
                str = get_bounds_vector (k, m, p);
                v_sum (b + str, column, -1, l);
            }
        }
        MPI_Barrier (MPI_COMM_WORLD);
    }
    //The Ux = b part of the solution
    MPI_Barrier (MPI_COMM_WORLD);
    if (l)
    {
        if (k % p == my_rank)
        {
            str = get_bounds(k, k, n, m, l, quan, re, my_rank, p);
            get_block(a, block, str, l, l);
            find_inv (k + 1, my_rank, p, n, m, l, quan, re,
                      a, block, inv, norrm, l);
            str2 = get_bounds_vector(k, m, p);
            get_vect (b, column, str2, str2 + l);
            mult (inv, l, l, 1, column, b + str2);
            get_vect(b, buf, str2, str2 + l);
            buf[l] = inv[l * l];
            MPI_Bcast (buf, l + 1, MPI_DOUBLE, my_rank, MPI_COMM_WORLD);
            if (inv[l * l] < 0)
            {
                return inv[l * l];
            }
        }
        else
        {
            MPI_Bcast (buf, l + 1, MPI_DOUBLE, k % p, MPI_COMM_WORLD);
            if (buf[l] < 0)
            {
                return buf[l];
            }
        }
        for (i = k - 1; i >= 0; i--)
        {
            if (my_rank == i % p)
            {
                str = get_bounds_vector (i, m, p);
                str2 = get_bounds (i, k, n, m, l, quan, re, my_rank, p);
                mult (a + str2, m, l, 1, buf, column);
                v_sum (b + str, column, -1, m);
            }
        }
    }
    MPI_Barrier (MPI_COMM_WORLD);
    for (j = k - 1; j >= 0; j--)
    {
        if (j % p == my_rank)
        {
            str = get_bounds(j, j, n, m, l, quan, re, my_rank, p);
            get_block(a, block, str, m, m);
            find_inv (j + 1, my_rank, p, n, m, l, quan, re,
                      a, block, inv, norrm, m);
            str2 = get_bounds_vector(j, m, p);
            get_vect (b, column, str2, str2 + m);
            mult (inv, m, m, 1, column, b + str2);
            get_vect(b, buf, str2, str2 + m);
            buf[m] = inv[m * m];
            MPI_Bcast (buf, m + 1, MPI_DOUBLE, my_rank, MPI_COMM_WORLD);
            if (inv[m * m] < 0)
            {
                return inv[m * m];
            }
        }
        else
        {
            MPI_Bcast (buf, m + 1, MPI_DOUBLE, j % p, MPI_COMM_WORLD);
            if (buf[m] < 0)
            {
                return buf[m];
            }
        }
        for (i = j - 1; i >= 0; i--)
        {
            if (my_rank == i % p)
            {
                str = get_bounds_vector (i, m, p);
                str2 = get_bounds (i, j, n, m, l, quan, re, my_rank, p);
                mult (a + str2, m, m, 1, buf, column);
                v_sum (b + str, column, -1, m);
            }
        }
        MPI_Barrier (MPI_COMM_WORLD);
    }
    return 0;
}
