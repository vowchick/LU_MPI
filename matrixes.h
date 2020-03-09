#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <mpi.h>
#define EPS 1e-16
#define max_print  10
#define N 60
struct args
{
    double *a = 0;                  //matrix
    double *x = 0;                  //answer
    double *b = 0;                  //rhs
    double *block = 0;
    int p = 0;                      // number of threads
    int n = 0;                      // matrix size
    int m = 0;                      //block size
    int id = 0;                     //thread id
    double cpu_time = 0;            
    double elapsed_time = 0;
    double add_time = 0;
    char *nname = 0;                //mutual filename
    int mist = 0;                   //in case of a mistake less than zero
    double *sum = 0;
};
void
gather_all (int j, int n, int m, int l, int quan, int re, int my_rank, int p,
                 int tot_quan, double *a, double *column, int *displs, int *recvcounts,
                 double *big_column, int m2);
void
send_vect (int j, int m, int p, double *b, double *buf, int my_rank);
void
mult_and_diff_v (int str, int col, int n, int m, int quan, int re, int my_rank, int p,
                 double *a, int left, int middle, int right, double *buf, double *column, double *b);
void
print_reg (double *a, int m, int n);
void
print_mistakes (char *name, int res);
double 
read_matrix (double *a, int n, int m, const char *name, int my_rank, int p);
double
init_matrix (double *a, int n, int m, int my_rank, int p);
void 
print_matrix (double *a, int n, int m, int my_rank, int p);
void
print_block_string (double *a, int n, int m_r, int m_c, int quan);
double 
f1(int i, int j, int n);
void
build_b (double *b, int n, double *a, int m, int my_rank, int p);
void
build_x (double *b, int n, int m, int my_rank, int p);
void
build_r (double *b, int n, int m, int my_rank, int p);
double 
norma_whole_matrix (double *a, int n, int m);
int
get_bounds1 (int i, int j, int n, int m);
int
solve (double *a, double *x, int n, int m, int my_rank, int p, double norrm, double *block);
int
Reverse (double *a, double *b, int n, int m, int my_rank,
        int p, double norrm, int re, int quan, double *block);
void 
prep_for_gather (int my_rank, int p, int quan, int *displs, int *recvcounts);
void find_inv (int j, int my_rank, int p, int n, int m, int l,
                         int quan, int re, double *a, 
                         double *block, double * inv, double norrm, int size);
int
get_bounds (int i/*строка*/, int j/*столбец*/, int n, int m, int l, 
            int quan, int re, int my_rank, int p);
int
get_bounds_vector (int i/*строка*/, int m, int p);
double
norm (double *y, double *z, int n);
double
init_forms (char *namea, double *a, int n, int m, int my_rank, int p);
double
Ax_b (double *a, double *x, double *b, int n, int m, 
      double *b_i, double *x_i, int p, int my_rank);
int
quantity (int j, int n, int m, int l, int quan, 
          int re, int my_rank, int p, int tot_quan, double *a, double *column, int m2);
void send_inv (int j, int my_rank, int p, int quan, int re, 
                      double *a, double *block, double *inv, double norrm, int m, int l, int n);
int 
calc_shift (int i, int p, int m, int m2, int *displs);
double 
AX_B_MPI (double *a, double *x,
          int n, int my_rank, int p, int m
          , double *b_i, double *x_i);
double
mult_str (double *str, double *x, double *b, int n, int m_c,
          int m_r, int p, int *displs);
int
LU (double *block, int b_size, double norrma);
double 
norma (double *matrix, int m_size);
int
LU_inv (double *block, double * inv, int b_size, double norrma);
void
delete_everything (double *a, double *b, double *c, double *d, double *e, double *f);
void
get_block (double *a, double *block, int start, int n1, int n2);
void
put_block (double *a, double *block, int start, int n1, int n2);
void
put_blocks_into_string (double *a, double *string, double *block, int str_num,
                        int m, int n, int quan, int re, int my_rank, int p);
void
put_block_into_string (double *str, double *block, int n, int m, int size1, int size2, int row_num);
void
put_inv (double *str, int start, double *block, int size);
int
get_bounds_row (int row_num, int row_col, int n, int m);
int
Inv_U (double *block, double * inv, int b_size, double norrma);
int
Inv_L (double *block, double * inv, int b_size, double norrma);
double 
norma_U (double *matrix, int m_size);
double 
norma_L (double *matrix, int m_size);
void 
mult (double *a, int m, int n, int k, double *b,double *c);
void
sum (double *a, double *b, int size); // видимо переписать
void
Reverse_Gauss (args *arg, double norrma);
void
v_sum (double *a, double *b, int coef, int size);
void
get_vect (double *b, double *vect, int i1, int i2);
void
put_vect (double *b, double *vect, int i1, int i2);
void 
sq_prod (double *a, double *b, double *c, int n);
void
sum2 (double *a, double *b, int coef, int size1, int size2);
void
put_zeros (double *a, int start, int size);
void
get_block_from_string (double *string, double *block, int i1, int j1, int i2, int j2, int m);
void
copy (double *dst, double *src, int size);
