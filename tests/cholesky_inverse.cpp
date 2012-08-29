#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include "mytimer.h"

inline float  __sqrt(const float  x) {return std::sqrt(x);}
inline double __sqrt(const double x) {return std::sqrt(x);}

#define MINEIG 0

extern "C" 
{
    void dpotrf_(char* U, int *N, double* A, int* lda, int* INFO);
    void dpotri_(char* U, int *N, double* A, int* lda, int* INFO);
    void dsyevd_(char*, char*, int*, double*, int *,double*, double*, int*, int*, int*, int*);
    
    void spotrf_(char* U, int *N, float* A, int* lda, int* INFO);
    void spotri_(char* U, int *N, float* A, int* lda, int* INFO);
    void ssyevd_(char*, char*, int*, float*, int *,float*, float*, int*, int*, int*, int*);

}

  template<int N, typename REAL>
REAL inverse_cholesky(const REAL in[N][N], REAL out[N][N])
{
  /* extract upper triangle from the [symmetric] input matrix */

  for (int j = 0; j < N; j++)
    for (int i = j; i < N; i++)
      out[j][i] = out[i][j] = in[j][i];

  /* Cholesky factorization, stored in the lower triangle */

  for (int k = 0; k < N; k++)
  {
    if (out[k][k] <= 0) return REAL(0.0);  /* matrix is not positive definite, return 0 */
    out[k][k] = __sqrt(out[k][k]);
    const REAL ainv = REAL(1.0)/out[k][k];
    for (int i = k+1; i < N; i++)
      out[i][k] *= ainv;
    for (int j = k+1; j < N; j++)
      for (int i = j; i < N; i++)
        out[i][j] -= out[i][k]*out[j][k];
  }

  for (int j = 0; j < N; j++)
    for (int i = j+1; i < N; i++)
      out[j][i] = 0.0;
  
  /* determinant */
  
  REAL det = 1.0;
  for (int i = 0; i < N; i++)
    det *= out[i][i];
  det *= det;
  assert(det > 0.0);

  /* invert lower triangular matrix */

  for (int k = 0; k < N; k++)
    out[k][k] = REAL(1.0)/out[k][k];

  for (int i = 1; i < N; i++)
    for (int j = 0; j < i; j++)
    {
      REAL sum = 0.0;
      for (int k = j; k < i; k++)
        sum += out[i][k] * out[k][j];
      out[i][j] = -out[i][i]*sum;
    }

  /* compute inverse by multiplying inverse of L and its transpose */

  for (int j = 0; j < N; j++)
    for (int i = j; i < N; i++)
    {
      REAL sum = 0;
      for (int k = i; k < N; k++)
        sum += out[k][j]*out[k][i];
      out[j][i] = sum;
    }

  for (int j = 0; j < N; j++)
    for (int i = j+1; i < N; i++)
      out[i][j] = out[j][i];

  return det;
}
  
  template<int N>
bool my_inverse(const double in[N][N], double out[N][N], double &det)
{
  int errorHandler;
  int     n = N;
  int lwork = N*N;
  char chU[] = "L";

  for (int j = 0; j < N; j++)
    for (int i = 0; i < N; i++)
      out[j][i] = in[j][i];

  dpotrf_(chU, &n, &out[0][0], &n, &errorHandler);
  assert(errorHandler >= 0);

  if (errorHandler > 0)
    return true;

  det = 1.0;
  for (int i = 0; i < N; i++)
    det *= out[i][i];
  det *= det;
  assert(det > 0.0);


  dpotri_(chU, &n, &out[0][0], &n, &errorHandler);
  assert(0 == errorHandler);
  for (int i = 0; i < N; i++)
    for (int j = i; j < N; j++)
      out[i][j] = out[j][i];

  return false;
}

  template<int N>
bool my_inverse(const float in[N][N], float out[N][N], float &det)
{
  int errorHandler;
  int     n = N;
  int lwork = N*N;
  char chU[] = "U";

  for (int j = 0; j < N; j++)
    for (int i = 0; i < N; i++)
      out[j][i] = in[j][i];

  spotrf_(chU, &n, &out[0][0], &n, &errorHandler);
  assert(errorHandler >= 0);

  if (errorHandler > 0)
    return true;

  det = 1.0;
  for (int i = 0; i < N; i++)
    det *= out[i][i];
  det *= det;
  assert(det > 0.0);


  spotri_(chU, &n, &out[0][0], &n, &errorHandler);
  assert(0 == errorHandler);
  for (int i = 0; i < N; i++)
    for (int j = i; j < N; j++)
      out[i][j] = out[j][i];

  return false;
}

  template<int N>
bool my_inverse_gold(const double in[N][N], double out[N][N], double &det)
{
  // pseudoinverse for real symmetric matrix
  const int N2 = N*N;
  double evalues[N];
  double evectors[N][N];

  for (int j = 0; j < N; j++)
    for (int i = 0; i < N; i++)
      evectors[j][i] = in[j][i];

  char jobz('V'), uplo('U');
  int   n = N;
  int lda = N;
  const int   LWORK = 1+6*N+2*N*N;
  const int  LIWORK = 3+5*N;

  double work[ LWORK];
  int   iwork[LIWORK];
  int   lwork =  LWORK;
  int  liwork = LIWORK;
  int info;

  dsyevd_(&jobz, &uplo, &n, &evectors[0][0], &lda, &evalues[0], &work[0],
      &lwork, &iwork[0], &liwork, &info);

  for (int i = 0; i < N; i++)
    for (int j = i; j < N; j++)
    {
      double recompinv = 0.0;
      for (int k = 0; k < N; k++)
        if (evalues[k] > MINEIG)
          recompinv += evectors[k][i]*evectors[k][j] / evalues[k];
      out[j][i] = out[i][j] = recompinv;
    }

  det = 1.0;

  bool filter = false;

  for (int i = 0; i < N; i++)
    if (evalues[i] < MINEIG)
    {
      filter  = true;
      det    *= MINEIG;
    }
    else
      det *= evalues[i];

  return filter;  
}

  template<int N>
bool my_inverse_gold(const float in[N][N], float out[N][N], float &det)
{
  // pseudoinverse for real symmetric matrix
  const int N2 = N*N;
  float evalues[N];
  float evectors[N][N];

  for (int j = 0; j < N; j++)
    for (int i = 0; i < N; i++)
      evectors[j][i] = in[j][i];

  char jobz('V'), uplo('U');
  int   n = N;
  int lda = N;
  const int   LWORK = 1+6*N+2*N*N;
  const int  LIWORK = 3+5*N;

  float work[ LWORK];
  int   iwork[LIWORK];
  int   lwork =  LWORK;
  int  liwork = LIWORK;
  int info;

  ssyevd_(&jobz, &uplo, &n, &evectors[0][0], &lda, &evalues[0], &work[0],
      &lwork, &iwork[0], &liwork, &info);

  for (int i = 0; i < N; i++)
    for (int j = i; j < N; j++)
    {
      float recompinv = 0.0;
      for (int k = 0; k < N; k++)
        if (evalues[k] > MINEIG)
          recompinv += evectors[k][i]*evectors[k][j] / evalues[k];
      out[j][i] = out[i][j] = recompinv;
    }

  det = 1.0;

  bool filter = false;

  for (int i = 0; i < N; i++)
    if (evalues[i] < MINEIG)
    {
      filter  = true;
      det    *= MINEIG;
    }
    else
      det *= evalues[i];

  return filter;  
}

template<int N, typename REAL>
void test(const int nrep)
{
  REAL mat[N][N] __attribute__((aligned(64)));
  REAL in [N][N] __attribute__((aligned(64)));
  REAL out[N][N] __attribute__((aligned(64)));
  REAL fdum;

  for (int j = 0; j < N; j++)
    for (int i = j; i < N; i++)
      mat[i][j] = mat[j][i] = drand48();

  for (int j = 0; j < N; j++)
    for (int i = 0; i < N; i++)
    {
      REAL res = 0.0;
      for (int k = 0; k < N; k++)
        res += mat[j][k] * mat[k][i];
      in[j][i] = res;
    }
  for (int j = 0; j < N; j++)
    for (int i = j; i < N; i++)
      assert(in[i][j] == in[j][i]);



  double t0 = get_wtime();
#pragma omp parallel for firstprivate(in) private(out)
  for (int r = 0; r < nrep; r++)
    assert(inverse_cholesky(in, out) > 0);
  const double dt_custom = get_wtime() - t0;
  fprintf(stderr, " custom inverse= %g sec \n", dt_custom);

  t0 = get_wtime();
#pragma omp parallel for firstprivate(in) private(out, fdum)
  for (int r = 0; r < nrep; r++)
    assert(my_inverse(in, out, fdum) == false);
  const double dt_lapack = get_wtime() - t0;
  fprintf(stderr, " LAPACK inverse= %g sec (ratio= %g)\n", dt_lapack, dt_lapack/dt_custom);

#if 0   /* does not work from MKL */
  t0 = get_wtime();
  for (int r = 0; r < nrep; r++)
    assert(my_inverse_gold(in, out, fdum) == false);
  fprintf(stderr, " DSYEVD inverse= %g sec \n", get_wtime() - t0);
#endif
};

template<typename REAL>
void do_test(const int n, const int nrep)
{
  switch(n)
  { 
    case  4:  test< 4, REAL>(nrep); break;
    case  5:  test< 5, REAL>(nrep); break;
    case  6:  test< 6, REAL>(nrep); break;
    case  7:  test< 7, REAL>(nrep); break;
    case  8:  test< 8, REAL>(nrep); break;
    case  9:  test< 9, REAL>(nrep); break;
    case 10:  test<10, REAL>(nrep); break;
    case 11:  test<11, REAL>(nrep); break;
    case 12:  test<12, REAL>(nrep); break;
    case 13:  test<13, REAL>(nrep); break;
    case 14:  test<14, REAL>(nrep); break;
    case 15:  test<15, REAL>(nrep); break;
    case 16:  test<16, REAL>(nrep); break;
    case 17:  test<17, REAL>(nrep); break;
    case 18:  test<18, REAL>(nrep); break;
    case 19:  test<19, REAL>(nrep); break;
    case 20:  test<20, REAL>(nrep); break;
    case 21:  test<21, REAL>(nrep); break;
    case 22:  test<22, REAL>(nrep); break;
    case 23:  test<23, REAL>(nrep); break;
    case 24:  test<24, REAL>(nrep); break;
    case 25:  test<25, REAL>(nrep); break;
    case 26:  test<26, REAL>(nrep); break;
    case 27:  test<27, REAL>(nrep); break;
    case 28:  test<28, REAL>(nrep); break;
    case 29:  test<29, REAL>(nrep); break;
    case 30:  test<30, REAL>(nrep); break;
    case 31:  test<31, REAL>(nrep); break;
    case 32:  test<32, REAL>(nrep); break;
    default: assert(0);
  }
}





int main(int argc, char *argv[])
{
  const int nrep =  argc > 1 ? atoi(argv[1]) : 100;
  fprintf(stderr, " Nrepeat= %d \n", nrep);
  const bool fp32 = argc > 2;
  fprintf(stderr, " Using %s precision \n", fp32 ? "SINGLE" : "DOUBLE" );

  /* warm-up phase */
  if (fp32)  do_test< float>(8, nrep);
  else       do_test<double>(8, nrep);



  for (int n = 4; n <= 32; n++)
  {
    fprintf(stderr, "Matrix size= %d x %d \n", n,n);
    if (fp32)  do_test< float>(n, nrep);
    else       do_test<double>(n, nrep);
  }

  const int N = 64;
  if (fp32) test<N, float>(nrep);
  else      test<N, double>(nrep);

  return 0;
}
