#include "matrixes.h"
double
norm (double * y, double *z, int n)
{
  double max = 0;
  int i;
  for (i = 0; i < n - 4; i += 4)
    {
      if (fabs (y[i] - z[i]) > max)
        max = fabs (y[i] - z[i]);
      if (fabs (y[i + 1] - z[i + 1]) > max)
        max = fabs (y[i + 1] - z[i + 1]);
      if (fabs (y[i + 2] - z[i + 2]) > max)
        max = fabs (y[i + 2] - z[i + 2]);
      if (fabs (y[i + 3] - z[i + 3]) > max)
        max = fabs (y[i + 3] - z[i + 3]);
    }
  for (; i < n; i++)
    if (fabs (y[i] - z[i]) > max)
        max = fabs (y[i] - z[i]);
  return max;
}
void
mult_and_diff_v (int str, int col, int n, int m, int quan, int re, int my_rank, int p,
                 double *a, int left, int middle, int right, double *buf, double *column, double *b)
{
    int str1 = 0, str2 = 0, l = n % m;
    str2 = get_bounds (str, col, n, m, l, quan, re, my_rank, p);
    mult (a + str2, left, middle, right, buf, column);
    str1 = get_bounds_vector (str, m, p);
    v_sum (b + str1, column, -1, left);
}
void
send_vect (int j, int m, int p, double *b, double *buf, int my_rank)
{
    int str = 0;
    str = get_bounds_vector (j, m, p);
    buf = b + str;
    MPI_Bcast (buf, m, MPI_DOUBLE, my_rank, MPI_COMM_WORLD);
}
int
quantity (int j, int n, int m, int l, int quan, 
          int re, int my_rank, int p, int tot_quan, double *a, double *column, int m2)
{
  int loc_q = 0, loc_quan = 0, str = 0, i = 0;  
  for (i = 0; i < j; i++)
    {
      if (my_rank == i % p)
        {
          str = get_bounds (i, j, n, m, l, quan, re, my_rank, p);
          if (loc_quan < tot_quan - l * m2)
            get_block (a, column + loc_q * m * m2, str, m, m2);
          else
            get_block (a, column + loc_q * m * m2, str, l, m2);
          loc_q++;
          if (loc_quan < tot_quan - l * m2)
            {
              loc_quan += m * m2;
            }
          else
            {
              loc_quan += l * m2;
            }
        }
    }
  return loc_quan;
}
int 
calc_shift (int i, int p, int m, int m2, int *displs)
{
  int shift = 0;
  shift = i / p;
  shift *= (m * m2); 
  shift += displs[i % p];
  return shift;
}
void send_inv (int j, int my_rank, int p, int quan, int re, 
                      double *a, double *block, double *inv, double norrm, int m, int l, int n)
{
  if (my_rank == (j - 1) % p)
    {
      find_inv (j, my_rank, p, n, m, l, quan, re, a,
                block, inv, norrm, m);
      MPI_Bcast (inv, m * m + 1, MPI_DOUBLE, my_rank, MPI_COMM_WORLD);
    }
  else
    {
      MPI_Bcast (inv, m * m + 1, MPI_DOUBLE, (j - 1) % p, MPI_COMM_WORLD);
    }
}
double 
AX_B_MPI (double *a, double *x,
          int n, int my_rank, int p, int m,
          double *b_i, double *x_i)
{
  int i, k = n / m, l = n % m;
  int quan = k / p, loc_q = 0;
  double max, max1 = 0, glob_max;
  max = 0;
  int re = k%p;
  int *recvcounts = new int[p];
  int *displs = new int [p];
  quan *= m;
  if (my_rank < re)
    {
        quan += m;
    }
  if(my_rank == re && l)
    {
        quan += l;
    }
  prep_for_gather (my_rank, p, quan, displs, recvcounts);
  MPI_Allgatherv (x_i, quan, MPI_DOUBLE, x, recvcounts, displs,
                  MPI_DOUBLE, MPI_COMM_WORLD);
  //printf ("my_rank = %d %d %d\n",my_rank, displs[my_rank], recvcounts[my_rank]);
  for (i = 0; i < k; i++)
    {
      if (my_rank == i % p)
        {
          max1 = mult_str (a + loc_q * n * m, x, b_i + loc_q * m, n, m, m, p, displs);
          if (max < max1)
            max = max1;
          loc_q++;
        }
    }
  if (my_rank == k % p)
    {
      max1 = mult_str (a + loc_q * n * m, x, b_i + loc_q * m, n, m, l, p, displs);
      if (max < max1)
            max = max1;
      loc_q++;
    }
  MPI_Allreduce (&max, &glob_max, 1, MPI_DOUBLE, 
                MPI_MAX, MPI_COMM_WORLD);
  delete []recvcounts;
  delete []displs;
  return glob_max;
}
void 
prep_for_gather (int my_rank, int p, int quan, int *displs, int *recvcounts)
{
  int i = 0, sen = 0, recc = 0, disp = 0, j = 0, rec = 0;
  int tag = 0;
  MPI_Status status;
  for (i = 0; i < p; i++)
    {
      if (my_rank < i)
        {
          sen = quan;
          MPI_Send (&sen, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
        }
      else if (my_rank == i)
        {
          disp = 0;
          for (j = 0; j < i; j++)
            {
              MPI_Recv (&recc, 1, MPI_INT, j, tag, MPI_COMM_WORLD, &status);
              disp += recc;
            }
        }
    }
  rec = quan;
  MPI_Allgather (&rec, 1, MPI_INT, recvcounts, 1, 
                MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather (&disp, 1, MPI_INT, displs, 1, 
                MPI_INT, MPI_COMM_WORLD);
}
int
get_bounds_row (int row_num, int row_col, int n, int m)
{
    int size1 = m, k = n / m, l = n % m;
    if (row_num == k)
    {
        size1 = l;
    }
    return row_col * size1 * m;
}
void find_inv (int j, int my_rank, int p, int n, int m, int l,
                         int quan, int re, double *a, 
                        double *block, double * inv, double norrm, int size)
{
  int str = 0, res = 0;
  str = get_bounds (j - 1, j - 1, n, m, l, quan,
              re, my_rank, p);
  //printf ("str1 = %d my_rank = %d \n", str, my_rank);
  get_block (a, block, str, size, size);
  res = LU_inv (block, inv, size, norrm);
  inv[size * size] = res;
}
double
mult_str (double *str, double *x, double *b, int n, int m_c,
          int m_r, int p, int *displs)
{
  int i, j, k = n / m_c, l = n % m_c, q, acc, indx = 0, indx2 = 0;
  double max = 0;
  double sum = 0;
  for (i = 0; i < m_r; i++)
    {
      sum = 0;
      for (j = 0; j < k; j++)
        {
          for (q = 0; q < m_c; q++)
            {
              acc = q + indx * m_c + displs[indx2];
              sum += str[q + j * m_c * m_r + i * m_c] * x[acc];
              //printf ("accdr = %d %.2f\n", q + j * m_c * m_r + i * n, str[q + j * m_c * m_r + i * n]);
            }
          if (indx2 < p - 1)
            {
              indx2++;
            }
          else
            {
              indx2 = 0;
              indx++;
            }
        }
      for (q = 0; q < l; q++)
        {
          acc = q + indx * m_c + displs[indx2];
          sum += str[q + k * m_c * m_r + i * l] * x[acc];
          //printf ("accdr = %d %.2f\n", q + k * m_c * m_r + i * n, str[q + k * m_c * m_r + i * n]);
        }
      //printf ("%d %.2f\n",i, b[i]);
      b[i] -= sum;
      b[i] *= (-1);
      if (max < fabs(b[i]))
        max = fabs (b[i]);
      sum = 0;
      indx = 0;
      indx2 = 0;
    }
  return max;
}
/*double
Ax_b (double *a, double *x, double *b, int n, int m, 
      double *b_i, double *x_i, int p, int my_rank) //переписать
{
  double sum = 0, max = 0;
  int tag = 0;
  MPI_Status status;
  int i, j, q, s, k = n / m, l = n - k * m;
  int k_p = k / p, k_pp = k % p;
  for (i = 0; i < k_p; i++)
    {
      for (s = 0; s < p; s++)
        {
          if (i * p + s == k)
            {
              for (j = 0; j < p; j++)
                MPI_Sendrecv (x_i, l, MPI_DOUBLE, j, tag, x + k * m,
                              l, MPI_DOUBLE, s, tag, MPI_COMM_WORLD, &status);
            }
          else
            {
              for (j = 0; j < p; j++)
                MPI_Sendrecv (x_i, m, MPI_DOUBLE, j, tag, x + i * p * m + s * m,
                              m, MPI_DOUBLE, s, tag, MPI_COMM_WORLD, &status);
            }
        }
    }
  for (i = 0; i < k_pp; i++)
    {
      if (
    }
  if (l != 0)
    {
      for (j = 0; j < p; j++)
        {

        }
    }
  for (q = 0; q < k; q++)
    {
      for (s = 0; s < m; s++)
        {
          for (j = 0; j < k; j++)
            {
              for (i = 0; i < m; i++)
                {
                  sum += a[i + j * m * m + s * m + q * m * n] * x[i + j * m];
                }
            }
          for (i = 0; i < l; i++)
                {
                  sum += a[i + k * m * m + s * l + q * m * n] * x[i + k * m];
                }
          sum -= b[q * m + s];
          if (fabs(sum) > max)
            max = fabs (sum);
          sum = 0;
        }
    }
    for (s = 0; s < l; s++)
      {
        for (j = 0; j < k; j++)
          {
            for (i = 0; i < m; i++)
              {
                sum += a[i + j * m * l + s * m + k * m * n] * x[i + j * m];
              }
          }
        for (i = 0; i < l; i++)
          {
            sum += a[i + k * m * l + s * l + k * m * n] * x[i + k * m];
          }
        sum -= b[k * m + s];
        if (fabs (sum) > max)
          max = fabs (sum);
        sum = 0;
      }
  return max;
}*/
int
get_bounds1 (int i, int j, int n, int m) // или переписать или добавить еще одну функцию
{
  int k = n / m, l = n - k * m;
  if (l != 0 && i == k)
    return n * m * i + j * m * l;
  return n * m * i + j * m * m;  
}
int
get_bounds (int i/*строка*/, int j/*столбец*/, int n, int m, int l, 
            int quan, int re, int my_rank, int p)
{
  int str = i / p;
  int block_st_i = str * (n) * m;
  if (str == (quan - 1) && my_rank == re && l)
    block_st_i += j * m * l;
  else
    block_st_i += j * m * m;
  return block_st_i;
}
int
get_bounds_vector (int i/*строка*/, int m, int p)
{
    int str = i / p;
    return str * m;
}
double 
norma_whole_matrix (double *a, int n, int m)  //не нужно вроде
{
  double sum = 0, max = 0;
  int i, j, q, s, k = n / m, l = n - k * m;
  for (q = 0; q < k; q++)
    {
      for (s = 0; s < m; s++)
        {
          for (j = 0; j < k; j++)
            {
              for (i = 0; i < m; i++)
                {
                  sum += fabs (a[i + j * m * m + s * m + q * m * n]);
                }
            }
          for (i = 0; i < l; i++)
                {
                  sum += fabs (a[i + k * m * m + s * l + q * m * n]);
                }
          if (fabs(sum) > max)
            max = fabs (sum);
          sum = 0;
        }
    }
    for (s = 0; s < l; s++)
      {
        for (j = 0; j < k; j++)
          {
            for (i = 0; i < m; i++)
              {
                sum += fabs (a[i + j * m * l + s * m + k * m * n]);
              }
          }
        for (i = 0; i < l; i++)
          {
            sum += fabs (a[i + k * m * l + s * l + k * m * n]);
          }
        if (fabs (sum) > max)
          max = sum;
        sum = 0;
      }
  return max;
}
int
LU (double *block, int b_size, double norrma)
{
  double norm = norma (block, b_size);
  if (norm < EPS * norrma)
    {
      return -1;
    }
  if (fabs (block[0]) < EPS * norm)
    return -1;
  int k, j;
  double rev = 1./block[0];
  for (k = 1; k < b_size - 4; k += 4)
    {
      block[k] *= rev;
      block[k + 1] *= rev;
      block[k + 2] *= rev;
      block[k + 3] *= rev;
    }
  for (; k < b_size; k++)
      block[k] *= rev;
  for (int i = 1; i < b_size; i++)
    {
      if (fabs (block[i * b_size + i]) < EPS * norm)
         return -1;
      for (k = 1; k <= i; k++)
        {
          double sum = 0;
          for (j = 0; j < k - 8; j+=8)
            {
              sum += block[i * b_size + j] * block[j * b_size + k];
              sum += block[i * b_size + j + 1] * block[(j + 1) * b_size + k];
              sum += block[i * b_size + j + 2] * block[(j + 2) * b_size + k];
              sum += block[i * b_size + j + 3] * block[(j + 3) * b_size + k];
              sum += block[i * b_size + j + 4] * block[(j + 4) * b_size + k];
              sum += block[i * b_size + j + 5] * block[(j + 5) * b_size + k];
              sum += block[i * b_size + j + 6] * block[(j + 6) * b_size + k];
              sum += block[i * b_size + j + 7] * block[(j + 7) * b_size + k];
            }
          for (; j < k; j++)
            sum += block[i * b_size + j] * block[j * b_size + k];
          block[i * b_size + k] -= sum;
        }
      for (k = b_size - 1; k > i; k--)
        {
          double sum = 0;
          for (j = 0; j < i - 8; j+=8)
            {
              sum += block[i * b_size + j] * block[j * b_size + k];
              sum += block[i * b_size + j + 1] * block[(j + 1) * b_size + k];
              sum += block[i * b_size + j + 2] * block[(j + 2) * b_size + k];
              sum += block[i * b_size + j + 3] * block[(j + 3) * b_size + k];
              sum += block[i * b_size + j + 4] * block[(j + 4) * b_size + k];
              sum += block[i * b_size + j + 5] * block[(j + 5) * b_size + k];
              sum += block[i * b_size + j + 6] * block[(j + 6) * b_size + k];
              sum += block[i * b_size + j + 7] * block[(j + 7) * b_size + k];
            }
          for (; j < i; j++)
            sum += block[i * b_size + j] * block[j * b_size + k];
          block[i * b_size + k] -= sum;
          if (fabs (block[i * b_size + i]) < EPS * norm)
             return -1;
          block[i * b_size + k] /= block[i * b_size + i]; 
        }
    }
  return 1;
}
double 
norma (double *matrix, int m_size)
{
  double max = 0;
  int j;
  for (int i = 0; i < m_size; i++)
    {
      double sum = 0;
      for (j = 0; j < m_size - 4; j += 4)
        {
          sum += fabs (matrix[i * m_size + j]);
          sum += fabs (matrix[i * m_size + j + 1]);
          sum += fabs (matrix[i * m_size + j + 2]);
          sum += fabs (matrix[i * m_size + j + 3]);
        }
      for (; j < m_size; j++)
        sum += fabs (matrix[i * m_size + j]);
      if (sum > max)
        max = sum;
    }
  return max;
}
int
LU_inv (double *block, double * inv, int b_size, double norrma)
{
  double norm = norma (block, b_size);
  if (norm < norrma * EPS)
    {
      return -1;
    }
  int res = LU (block, b_size, norrma), k;
  if (res < 0)
    return res;
  for (int i = 0; i < b_size; i++)
    {
      double sum;
      for (int j = 0; j < b_size; j++)
        {
          if (fabs (block[j * b_size + j]) < EPS * norm)
            return -1;
          sum = 0;
          for (k = 0; k < j - 8; k += 8)
            {
              sum += inv[k * b_size + i] * block[j * b_size + k];
              sum += inv[(k + 1) * b_size + i] * block[j * b_size + k + 1];
              sum += inv[(k + 2) * b_size + i] * block[j * b_size + k + 2];
              sum += inv[(k + 3) * b_size + i] * block[j * b_size + k + 3];
              sum += inv[(k + 4) * b_size + i] * block[j * b_size + k + 4];
              sum += inv[(k + 5) * b_size + i] * block[j * b_size + k + 5];
              sum += inv[(k + 6) * b_size + i] * block[j * b_size + k + 6];
              sum += inv[(k + 7) * b_size + i] * block[j * b_size + k + 7];
            }
          for (;k < j; k++)
            sum += inv[k * b_size + i] * block[j * b_size + k];
          inv[j * b_size + i] = ((i == j ? 1 : 0) - sum) / block[j * b_size + j];
        }
      for (int j = b_size - 1; j >= 0; j--)
        {
          sum = 0;
          for (int k = b_size - 1; k > j; k--)
            sum += inv[k * b_size + i] * block[j * b_size + k];
          inv[j * b_size + i] -= sum;
        }
    }
  return 1;
}
void
get_block (double *a, double *block, int start, int n1/*строки*/, int n2/*столбцы*/)
{
  int i, j;
  for (i = 0; i < n1; i++)
    {
      for (j = 0; j < n2 - 4; j += 4)
        {
            block[i * n2 + j] = a[start + i * n2 + j];
            block[i * n2 + j + 1] = a[start + i * n2 + j + 1];
            block[i * n2 + j + 2] = a[start + i * n2 + j + 2];
            block[i * n2 + j + 3] = a[start + i * n2 + j + 3];
        }
      for  (;j < n2; j++)
        {
            block[i * n2 + j] = a[start + i * n2 + j];
        }
    }
}
void
put_block (double *a, double *block, int start, int n1/*строки*/, int n2/*столбцы*/)
{
  int i, j;
  for (i = 0; i < n1; i++)
    {
      for (j = 0; j < n2; j++)
        {
            a[start + i * n2 + j] = block[i * n2 + j];
        }
    }
}
int
Inv_L (double *block, double * inv, int b_size, double norrma)
{
  double norm = norma_L (block, b_size);
  if (norm < norrma * EPS)
  {
    return -1;
  }
  int k;
  for (int i = 0; i < b_size; i++)
    {
      double sum;
      for (int j = 0; j < b_size; j++)
        {
          if (fabs (block[j * b_size + j]) < EPS * norm)
            return -1;
          sum = 0;
          for (k = 0; k < j - 8; k+=8)
            {
              sum += inv[k * b_size + i] * block[j * b_size + k];
              sum += inv[(k + 1) * b_size + i] * block[j * b_size + k + 1];
              sum += inv[(k + 2) * b_size + i] * block[j * b_size + k + 2];
              sum += inv[(k + 3) * b_size + i] * block[j * b_size + k + 3];
              sum += inv[(k + 4) * b_size + i] * block[j * b_size + k + 4];
              sum += inv[(k + 5) * b_size + i] * block[j * b_size + k + 5];
              sum += inv[(k + 6) * b_size + i] * block[j * b_size + k + 6];
              sum += inv[(k + 7) * b_size + i] * block[j * b_size + k + 7];
            }
          for (; k < j; k++)
            sum += inv[k * b_size + i] * block[j * b_size + k];
          if (fabs (block[j * b_size + j]) < EPS * norm)
            return -1;
          inv[j * b_size + i] = ((i == j ? 1 : 0) - sum) / block[j * b_size + j];
        }
    }
  return 1;
}
double 
norma_U (double *matrix, int m_size)
{
  double max = 0;
  int j;
  for (int i = 0; i < m_size; i++)
    {
      double sum = 1;
      for (j = i + 1; j < m_size - 4; j+=4)
        {
          sum += fabs (matrix[i * m_size + j]);
          sum += fabs (matrix[i * m_size + j + 1]);
          sum += fabs (matrix[i * m_size + j + 2]);
          sum += fabs (matrix[i * m_size + j + 3]);
        }
      for (; j < m_size; j++)
        sum += fabs (matrix[i * m_size + j]);
      if (sum > max)
        max = sum;
    }
  return max;
}
double 
norma_L (double *matrix, int m_size)
{
  double max = 0;
  int j = 0;
  for (int i = 0; i < m_size; i++)
    {
      double sum = 0;
      for (j = 0; j <= i - 4; j+=4)
        {
          sum += fabs (matrix[i * m_size + j]);
          sum += fabs (matrix[i * m_size + j + 1]);
          sum += fabs (matrix[i * m_size + j + 2]);
          sum += fabs (matrix[i * m_size + j + 3]);
        }
      for (; j <= i; j++)
        sum += fabs (matrix[i * m_size + j]);
      if (sum > max)
        max = sum;
    }
  return max;
}
int
Inv_U (double *block, double * inv, int b_size, double norrma)
{
  double norm = norma_U (block, b_size);
  if (norm < norrma * EPS)
  {
    return -1;
  }
  int k;
  for (int i = 0; i < b_size; i++)
    {
      double sum;
      for (int j = b_size - 1; j >= 0; j--)
        {
          sum = 0;
          for (k = b_size - 1; k > j + 8; k -= 8)
            {
              sum += inv[k * b_size + i] * block[j * b_size + k];
              sum += inv[(k - 1) * b_size + i] * block[j * b_size + k - 1];
              sum += inv[(k - 2) * b_size + i] * block[j * b_size + k - 2];
              sum += inv[(k - 3) * b_size + i] * block[j * b_size + k - 3];
              sum += inv[(k - 4) * b_size + i] * block[j * b_size + k - 4];
              sum += inv[(k - 5) * b_size + i] * block[j * b_size + k - 5];
              sum += inv[(k - 6) * b_size + i] * block[j * b_size + k - 6];
              sum += inv[(k - 7) * b_size + i] * block[j * b_size + k - 7];
            }
          for (; k > j; k--)
            sum += inv[k * b_size + i] * block[j * b_size + k];
          inv[j * b_size + i] = ((i == j ? 1 : 0) - sum);
        }
    }
  return 1;
}
void 
mult (double *a, int m, int n, int k, double *b,double *c)
{
  double sum = 0;
  int z = n & 3;
  int j;
  for (int p = 0; p < m; p++)
    {
      for (int i = 0; i < k; i++)
        {
          for (j = 0; j < z; j++)
            sum += a[p * n + j] * b[j * k + i];
          for (j = z; j < n; j += 4)
            {
              sum += a[p * n + j] * b[j * k + i];
              sum += a[p * n + j + 1] * b[(j + 1) * k + i];
              sum += a[p * n + j + 2] * b[(j + 2) * k + i];
              sum += a[p * n + j + 3] * b[(j + 3) * k + i];
            }
          c[p * k + i] = sum;
          sum = 0;
        }
    }
}
void
sum2 (double *a, double *b, int coef, int size1, int size2)
{
    for (int i = 0; i < size1; i++)
    {
      for (int j = 0; j < size2; j++)
        a[i * size2 + j] += (b[i * size2 + j] * coef);
    }
}
void
gather_all (int j, int n, int m, int l, int quan, int re, int my_rank, int p,
                 int tot_quan, double *a, double *column, int *displs, int *recvcounts,
                 double *big_column, int m2)
{
  int loc_quan = 0;
  loc_quan = quantity (j, n, m, l, quan, re, my_rank, p,
                      tot_quan, a, column, m2);
  prep_for_gather (my_rank, p, loc_quan, displs, recvcounts);
  MPI_Allgatherv (column, loc_quan, MPI_DOUBLE, big_column, recvcounts,
                  displs, MPI_DOUBLE, MPI_COMM_WORLD);
}
void
sum (double *a, double *b, int size)
{
  int j;
  for (int i = 0; i < size; i++)
  {
    for (j = 0; j < size - 8; j+=8)
      {
        a[i * size + j] -= (b[i * size + j]);
        a[i * size + j + 1] -= (b[i * size + j + 1]);
        a[i * size + j + 2] -= (b[i * size + j + 2]);
        a[i * size + j + 3] -= (b[i * size + j + 3]);
        a[i * size + j + 4] -= (b[i * size + j + 4]);
        a[i * size + j + 5] -= (b[i * size + j + 5]);
        a[i * size + j + 6] -= (b[i * size + j + 6]);
        a[i * size + j + 7] -= (b[i * size + j + 7]);
      }
    for (; j < size; j++)
      a[i * size + j] -= (b[i * size + j]);
  }
}
void
v_sum (double *a, double *b, int coef, int size)
{
  for (int i = 0; i < size; i++)
    a[i] += (b[i] * coef);
}
void
get_vect (double *b, double *vect, int i1, int i2)
{
  for (int i = 0; i < i2 - i1; i++)
    vect[i] = b[i1 + i];
}
void
put_vect (double *b, double *vect, int i1, int i2)
{
  for (int i = 0; i < i2 - i1; i++)
    b[i1 + i] = vect[i];
}
void
sq_prod1 (double *a, double *b, double *c, int n)
{
  int bm, bi, nbm, nbi;
  int l, nl;
  int i, j, m;
  double *pa, *pb, *pc;
  double s00, s01, s02, s03;
  double s10, s11, s12, s13;
  double s20, s21, s22, s23;
  double s30, s31, s32, s33;
  for (bm = 0; bm < n; bm += N)
    {
      nbm = (bm + N <= n ? bm + N : n);
      for (bi = 0; bi < n; bi += N)
        {
          nbi = (bi + N <= n ? bi + N : n);
          for (m = bm, pc = c + bm; m < nbm; m++, pc++)
            for (i = bi; i < nbi; i++)
              pc[i * n] = 0.;
          for (l = 0; l < n; l += N)
            {
              nl = (l + N <= n ? l + N : n);
              for (m = bm, pc = c + bm; m < nbm;
                   m += 4, pc += 4)
              for (i = bi, pb = b + m; i < nbi; i += 4)
                {
                  pa = a + l + i * n;
                  s00 = s01 = s02 = s03 = 0.;
                  s10 = s11 = s12 = s13 = 0.;
                  s20 = s21 = s22 = s23 = 0.;
                  s30 = s31 = s32 = s33 = 0.;
                  for (j = l; j < nl; j++, pa++)
                    {
                      s00 += pa[0] * pb[j * n];
                      s01 += pa[0] * pb[j * n + 1];
                      s02 += pa[0] * pb[j * n + 2];
                      s03 += pa[0] * pb[j * n + 3];

                      s10 += pa[n] * pb[j * n];
                      s11 += pa[n] * pb[j * n + 1];
                      s12 += pa[n] * pb[j * n + 2];
                      s13 += pa[n] * pb[j * n + 3];

                      s20 += pa[2 * n] * pb[j * n];
                      s21 += pa[2 * n] * pb[j * n + 1];
                      s22 += pa[2 * n] * pb[j * n + 2];
                      s23 += pa[2 * n] * pb[j * n + 3];

                      s30 += pa[3 * n] * pb[j * n];
                      s31 += pa[3 * n] * pb[j * n + 1];
                      s32 += pa[3 * n] * pb[j * n + 2];
                      s33 += pa[3 * n] * pb[j * n + 3];
                    }
                  pc[i * n]             += s00;
                  pc[i * n + 1]         += s01;
                  pc[i * n + 2]         += s02;
                  pc[i * n + 3]         += s03;

                  pc[(i + 1) * n]       += s10;
                  pc[(i + 1) * n + 1]   += s11;
                  pc[(i + 1) * n + 2]   += s12;
                  pc[(i + 1) * n + 3]   += s13;

                  pc[(i + 2) * n]       += s20;
                  pc[(i + 2) * n + 1]   += s21;
                  pc[(i + 2) * n + 2]   += s22;
                  pc[(i + 2) * n + 3]   += s23;

                  pc[(i + 3) * n]       += s30;
                  pc[(i + 3) * n + 1]   += s31;
                  pc[(i + 3) * n + 2]   += s32;
                  pc[(i + 3) * n + 3]   += s33;
                }
            }
        }
    }
}
void
put_zeros (double *a, int start, int size)
{
  for (int i = 0; i < size; i++)
    {
      a[i + start] = 0.;
    }
}
void
put_blocks_into_string (double *a, double *string, double *block, int str_num,
                        int m, int n, int quan, int re, int my_rank, int p)
{
  int k = n / m, l = n % m;
  int i = 0, j = 0, str = 0, size1 = 0;
  if (str_num == k)
      size1 = l;
  else
      size1 = m;
  for (i = 0; i < k; i++)
    {
        str = get_bounds (str_num, i, n, m, l, quan, re, my_rank, p);
        get_block (a, block, str, size1, m);
        put_block_into_string (string, block, n, m, size1, m, i);
    }
  if (l)
    {
      put_block_into_string (string, block, n, m, size1, l, k);
    }
}
void
put_block_into_string (double *str, double *block, int n, int m, int size1, int size2, int row_num)
{
  int i = 0, j = 0;
  for (i = 0; i < size1; i++)
    {
      for (j = 0; j < size2; j++)
        {
          str[row_num * size1 * size2 + j + i * size2] = block[i * size2 + j];
        }
    }
}
void
put_inv (double *str, int start, double *block, int size)
{
    int i = 0,j = 0;
    for (i = 0; i < size; i++)
        {
            for (j = 0; j < size; j++)
                str[start + i * size + j] = block[i * size + j];
        }
    str[start + size * size] = block[size * size];
}
void
get_block_from_string (double *string, double *block, int i1, int j1, int i2, int j2, int m)
{
    int height = i2 - i1, length = j2 - j1, j;
    for (int i = 0; i < height; i++)
      {
        for (j = 0; j < length - 4; j += 4)
          {
            block[i * length + j] = string[(i1 + i) * m + (j1 + j)];
            block[i * length + j + 1] = string[(i1 + i) * m + (j1 + j + 1)];
            block[i * length + j + 2] = string[(i1 + i) * m + (j1 + j + 2)];
            block[i * length + j + 3] = string[(i1 + i) * m + (j1 + j + 3)];
          }
        for (; j < length; j++)
          block[i * length + j] = string[(i1 + i) * m + (j1 + j)];
      }
}
void
copy (double *a, double *b, int n)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            a[i * n + j] = b[i * n + j];
}
void
sq_prod (double *a, double *b, double *c, int n)
{
  int bm, bi, nbm, nbi;
  int l, nl;
  int i, j, m;
  double *pa, *pb, *pc;
  double s00, s01, s02;
  double s10, s11, s12;
  double s20, s21, s22;
  for (bm = 0; bm < n; bm += N)
    {
      nbm = (bm + N <= n ? bm + N : n);
      for (bi = 0; bi < n; bi += N)
        {
          nbi = (bi + N <= n ? bi + N : n);
          for (m = bm, pc = c + bm; m < nbm; m++, pc++)
            for (i = bi; i < nbi; i++)
              pc[i * n] = 0.;
          for (l = 0; l < n; l += N)
            {
              nl = (l + N <= n ? l + N : n);
              for (m = bm, pc = c + bm; m < nbm;
                   m += 3, pc += 3)
              for (i = bi, pb = b + m; i < nbi; i += 3)
                {
                  pa = a + l + i * n;
                  s00 = s01 = s02 = 0.;
                  s10 = s11 = s12 = 0.;
                  s20 = s21 = s22 = 0.;
                  for (j = l; j < nl; j++, pa++)
                    {
                      s00 += pa[0] * pb[j * n];
                      s01 += pa[0] * pb[j * n + 1];
                      s02 += pa[0] * pb[j * n + 2];

                      s10 += pa[n] * pb[j * n];
                      s11 += pa[n] * pb[j * n + 1];
                      s12 += pa[n] * pb[j * n + 2];

                      s20 += pa[2 * n] * pb[j * n];
                      s21 += pa[2 * n] * pb[j * n + 1];
                      s22 += pa[2 * n] * pb[j * n + 2];
                    }
                  pc[i * n]             += s00;
                  pc[i * n + 1]         += s01;
                  pc[i * n + 2]         += s02;

                  pc[(i + 1) * n]       += s10;
                  pc[(i + 1) * n + 1]   += s11;
                  pc[(i + 1) * n + 2]   += s12;

                  pc[(i + 2) * n]       += s20;
                  pc[(i + 2) * n + 1]   += s21;
                  pc[(i + 2) * n + 2]   += s22;

                }
            }
        }
    }
}
