#include <assert.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>

unsigned int m=0;
unsigned int n=0;
double **X;


/*
 * create a new n1 x n2 matrix which is allocated like
 * an n1*n2 array, but can be referenced as a 2-d array
 */

double ** new_matrix(unsigned int n1, unsigned int n2)
{
  int i;
  double **m;
  
  if(n1 == 0 || n2 == 0) return NULL;
  
  m = (double**) malloc(sizeof(double*) * n1);
  assert(m);
  m[0] = (double*) malloc(sizeof(double) * (n1*n2));
  assert(m[0]);
  
  for(i=1; i<n1; i++) m[i] = m[i-1] + n2;
  
  return m;
}


/*
 * create a double ** Matrix from a double * vector
 * should be freed with the free command, rather than
 * delete_matrix
 */

double ** new_matrix_bones(double *v, unsigned int n1, unsigned int n2)
{
  double **M;
  int i;
  M = (double **) malloc(sizeof(double*) * n1);
  M[0] = v;
  for(i=1; i<n1; i++) M[i] = M[i-1] + n2;
  return(M);
}

/*
 * copies vold to v
 * (assumes v has already been allcocated)
 */

void dupv(double *v, double* vold, unsigned int n)
{
  unsigned int i;
  for(i=0; i<n; i++) v[i] = vold[i];
}


/*
 * copy M2 to M1
 */

void dup_matrix(double** M1, double **M2, unsigned int n1, unsigned int n2)
{
  unsigned int i;
  if(n1 == 0 || n2 == 0) return;
  assert(M1 && M2);
  for(i=0; i<n1; i++) dupv(M1[i], M2[i], n2);
}

/*
 * create a new n1 x n2 matrix which is allocated like
 * an n1*n2 array, and copy the of n1 x n2 M into it.
 */

double ** new_dup_matrix(double** M, unsigned int n1, unsigned int n2)
{
  double **m;

  if(n1 <= 0 || n2 <= 0) {
    /* assert(M == NULL); */
    return NULL;
  }

  m = new_matrix(n1, n2);
  dup_matrix(m, M, n1, n2);
  return m;
}


/*
 * delete a matrix allocated as above
 */

void delete_matrix(double** m)
{
  if(m == NULL) return;
  assert(*m);
  free(*m);
  assert(m);
  free(m);
}


/*
 * sq:
 * 
 * calculate the square of x
 */

double sq(double x)
{
  return x*x;
}


/*
 * min_of_columns:
 *
 * fill s[n1] with the min of the columns of M (n1 x n2);
 */

void min_of_columns(double *s, double **M, unsigned int n1, unsigned int n2)
{
  unsigned int i,j;

  /* sanity checks */
  if(n1 <= 0 || n2 <= 0) {return;}
  assert(s && M);
  
  /* calculate sum of columns */
  for(i=0; i<n2; i++) {
    s[i] = M[0][i];
    for(j=1; j<n1; j++) if(M[j][i] < s[i]) s[i] = M[j][i];
  }
}


/*
 * new vector of integers of length n
 */

int *new_ivector(unsigned int n)
{
  int *iv;
  if(n == 0) return NULL;
  iv = (int*)  malloc(sizeof(int) * n);
  assert(iv);
  return iv;
}


/*
 * allocate and return an array containing
 * the integer seqence [from...to]
 */

int* iseq(double from, double to)
{
  unsigned int n,i;
  int by;
  int *s = NULL;
  
  if(from <= to) {
    n = (unsigned int) (to - from) + 1;
    by = 1;
  } else {
    assert(from > to);
    n = (unsigned int) (from - to) + 1;
    by = -1;
  }
  
  if(n == 0) return NULL;
  
  s = new_ivector(n);
  s[0] = from;
  for(i=1; i<n; i++) {
    s[i] = s[i-1] + by;
  }
  return s;
}


/*
 * duplicate the integer contents of iv of length n into the already
 * allocated vector iv_new, also of length n
 */

void dupiv(int *iv_new, int *iv, unsigned int n)
{
  unsigned int i;
  if(n > 0) assert(iv && iv_new);
  for(i=0; i<n; i++) iv_new[i] = iv[i];
}


/*
 * distance:
 *
 * C-side version of distance_R
 */

void distance(double **X1, const unsigned int n1, double **X2,
              const unsigned int n2, const unsigned int m,
              double **D)
{
  unsigned int i,j,k;

  /* for each row of X1 and X2 */
  for(i=0; i<n1; i++) {
    for(j=0; j<n2; j++) {

      /* sum the squared entries */
      D[i][j] = 0.0;
      for(k=0; k<m; k++) {
        D[i][j] += sq(X1[i][k] - X2[j][k]);
      }

    }
  }
}

/*
 * Returns the kth smallest value in the array arr[1..n].  The input
 * array will be rearranged to have this value in location arr[k] ,
 * with all smaller elements moved to arr[1..k-1] (in arbitrary order)
 * and all larger elements in arr[k+1..n] (also in arbitrary order).
 * (from Numerical Recipies in C)
 *
 * This Quickselect routine is based on the algorithm described in
 * "Numerical recipes in C", Second Edition, Cambridge University
 * Press, 1992, Section 8.5, ISBN 0-521-43108-5 This code by Nicolas
 * Devillard - 1998. Public domain.
 */

#define ELEM_SWAP(a,b) { register double t=(a);(a)=(b);(b)=t; }
#define IELEM_SWAP(a,b) { register int t=(a);(a)=(b);(b)=t; }

double quick_select_index(double arr[], int iarr[], int n, int k)
{
  int low, high ;
  int middle, ll, hh;

  low = 0 ; high = n-1 ;
  assert(k >= low && k <= high);
  for (;;) {
    if (high <= low) /* One element only */
      return arr[k] ;

    if (high == low + 1) {  /* Two elements only */
      if (arr[low] > arr[high]) {
        ELEM_SWAP(arr[low], arr[high]) ;
        if(iarr) IELEM_SWAP(iarr[low], iarr[high]) ;
      }
      return arr[k] ;
    }

    /* Find kth of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if (arr[middle] > arr[high]) {
      ELEM_SWAP(arr[middle], arr[high]) ;
      if(iarr) IELEM_SWAP(iarr[middle], iarr[high]) ;
    }
    if (arr[low] > arr[high]) {
      ELEM_SWAP(arr[low], arr[high]) ;
      if(iarr) IELEM_SWAP(iarr[low],iarr[high]) ;
    }
    if (arr[middle] > arr[low])  {
      ELEM_SWAP(arr[middle], arr[low]) ;
      if(iarr) IELEM_SWAP(iarr[middle],iarr[low]) ;
    }

    /* Swap low item (now in position middle) into position (low+1) */
    ELEM_SWAP(arr[middle], arr[low+1]) ;
    if(iarr) IELEM_SWAP(iarr[middle], iarr[low+1]) ;

    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
        do ll++; while (arr[low] > arr[ll]) ;
        do hh--; while (arr[hh]  > arr[low]) ;

        if (hh < ll)
        break;

        ELEM_SWAP(arr[ll], arr[hh]) ;
        if(iarr) IELEM_SWAP(iarr[ll], iarr[hh]) ;
    }

    /* Swap middle item (in position low) back into correct position */
    ELEM_SWAP(arr[low], arr[hh]) ;
    if(iarr) IELEM_SWAP(iarr[low], iarr[hh]) ;

    /* Re-set active partition */
    if (hh <= k)
        low = ll;
    if (hh >= k)
        high = hh - 1;
  }
}
/*
 * closest:
 *
 * returns the start indices into X which are closest (via min for nref
 * > 1) to the element(s) of Xref.  The first start of those indices are
 * the start closest, otherwise the indices unordered (unless sorted=true).
 * Even when sorted=true the indices close+1, ... are not sorted.
 */

int *closest(const unsigned int m, const unsigned int start,
  double **Xref, const unsigned int nref, const unsigned int n, double **X)
{
  int *oD;
  double **D;

  /* calculate distances to reference location(s), and so-order X & Z */
  D = new_matrix(nref, n);
  distance(Xref, nref, X, n, m, D);
  if(nref > 1) min_of_columns(*D, D, nref, n);

  /* partition based on "close"st */
  if(n > start) {
    oD = iseq(1, n);
    quick_select_index(*D, oD, n, start);
  } else oD = NULL;

  delete_matrix(D);
  return(oD);
}

/*
 * closest_R:
 *
 * R interface to closest indices, primarily for debugging purposes
 */

void closest_R(

  /* inputs */
  int *m_in, 
  int *start_in,
  double *Xref_in,
  int *nref_in,
  /* double *X_in, */

  /* outputs */
  int *oD_out)

{
  double /* **X, */ **Xref;
  int *oD;

  /* make matrix bones */
  /* X = new_matrix_bones(X_in, *n_in, *m_in); */
  assert(*start_in < n);
  assert(*m_in == m);
  Xref = new_matrix_bones(Xref_in, *nref_in, *m_in);

  oD = closest(m, *start_in, Xref, *nref_in, n, X);

  dupiv(oD_out, oD, *start_in);
  free(oD); 

  /* free(X); */
  free(Xref);
}


/* 
 * loadX_R:
 *
 * bulk copy X over to this side so that multiple 
 * closest calls on it can be fast via pointers
 */

void loadX_R(int *m_in, int *n_in, double *X_in)
{
  double **Xtemp;
  if(X) delete_matrix(X);
  Xtemp = new_matrix_bones(X_in, *n_in, *m_in);
  n = *n_in;
  m = *m_in;
  X = new_dup_matrix(Xtemp, n, m);
  free(Xtemp);
}


/* 
 * unloadX_R:
 *
 * clean up 
 */

void unloadX_R(void)
{
  if(X) delete_matrix(X);
  m = n = 0;
  X = NULL;
}
