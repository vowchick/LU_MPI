#include "matrixes.h"
//#include "time.h"
//in case of a mistake should return -10
int
solve (double *a, double *b, double *x, int n, int m, int my_rank, int p, double norrm)
{
    int i, j, s, k = n / m, l = n % m;
    int str = 0, shift = 0, str2 = 0;
    int *recvcounts = new int[2 * p];
    int *displs = recvcounts + p;
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
    double *block, *inv, *column, *big_column;
    block = new double [2 * m * m + 1 + tot_quan + n * m];
    inv = block + m * m;
    column = inv + m * m + 1;
    big_column = column + tot_quan;
    for (j = 1; j < k; j++)
    {
        /*
        for (s = j; j < n; s++)
          L[s][j − 1] = A[s][j − 1] ∗ inv (A[j − 1][j − 1]);
      */
        send_inv (j, my_rank, p, quan, re, a, block, inv, norrm, m, l, n);
        if (inv[m * m] < 0)
        {
            delete []block;
            delete []recvcounts;
            return inv[m * m];
        }
        for (s = j; s < k; s++)
        {
            if (my_rank == s % p)
            {
                str = get_bounds (s, j - 1, n, m, l, quan, re, my_rank, p);
                get_block (a, block, str, m, m);
                sq_prod (block, inv, a + str, m);
            }
        }
        if (l != 0 && (my_rank == (k % p)))
        {
            str = get_bounds (k, j - 1, n, m, l, quan, re, my_rank, p);
            get_block (a, block, str, l, m);
            mult (block, l, m, m, inv, a + str);
        }
        // теперь пересылка всем процессам столбца j от 0 до j
        gather_all (j, n, m, l, quan, re, my_rank, p, tot_quan, a, column,
                    displs, recvcounts, big_column, m);
        for (i = 0; i < j; i++)
        {
            shift = calc_shift (i, p, m, m, displs);
            for (s = i + 1; s < k; s++)
            {
                //A [ s ] [ j ] −= L [ s ] [ i ] ∗ A [ i ] [ j ]

                if (my_rank == s % p)
                {
                    str = get_bounds (s, i, n, m, l, quan, re, my_rank, p);
                    sq_prod (a + str, big_column + shift, block, m);
                    str2 = get_bounds (s, j, n, m, l, quan, re, my_rank, p);
                    sum (a + str2, block, m);
                }
            }
            if (l && (my_rank == k % p))
            {
                str = get_bounds (k, i, n, m, l, quan, re, my_rank, p);
                mult (a + str, l, m, m, big_column + shift, block);
                str = get_bounds (k, j, n, m, l, quan, re, my_rank, p);
                sum2 (a + str, block, -1, l, m);
            }
        }
    }
    if (l)
    {
        send_inv (k, my_rank, p, quan, re, a, block, inv, norrm, m, l, n);
        if (inv[m * m] < 0)
        {
            delete []block;
            delete []recvcounts;
            return inv[m * m];
        }
        if (my_rank == k % p)
        {
            str = get_bounds (k, k - 1, n, m, l, quan, re, my_rank, p);
            get_block (a, block, str, l, m);
            mult (block, l, m, m, inv, a + str);
        }
        gather_all (k, n, m, l, quan, re, my_rank, p, tot_quan, a, column,
                    displs, recvcounts, big_column, l);
        /*for (i = 0; i < j; i++)
          {
            for (int q = 0; q < m; q++)
              {
                for (int s = 0; s < l; s++)
                  printf ("%.2f ", big_column[i * m * l + s + q * l]);
                printf ("\n");
              }
          }*/
        for (i = 0; i < k; i++)
        {
            shift = calc_shift (i, p, m, l, displs);
            for (s = i + 1; s < k; s++)
            {
                //A [ s ] [ j ] −= L [ s ] [ i ] ∗ A [ i ] [ j ]
                
                if (my_rank == s % p)
                {
                    str = get_bounds (s, i, n, m, l, quan, re, my_rank, p);
                    mult (a + str, m, m, l, big_column + shift, block);
                    str = get_bounds (s, k, n, m, l, quan, re, my_rank, p);
                    sum2 (a + str, block, -1, m, l);
                }
            }
            if (my_rank == k % p)
            {
                str = get_bounds (k, i, n, m, l, quan, re, my_rank, p);
                mult (a + str, l, m, l, big_column + shift, block);
                str = get_bounds (k, k, n, m, l, quan, re, my_rank, p);
                sum (a + str, block, l);
            }
        }
    }
    MPI_Barrier (MPI_COMM_WORLD);
    print_matrix(a, n, m, my_rank, p);
    MPI_Barrier (MPI_COMM_WORLD);
    Reverse (a, b, x, n, m, my_rank, p, norrm, re, quan, block);
    delete []block;
    delete []recvcounts;
    return 0;
}
int
Reverse (double *a, double *b, double *x,
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
