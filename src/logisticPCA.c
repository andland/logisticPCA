#include <R.h>
#include <Rinternals.h>
#include <math.h>


static inline double logistic(double x)
{
  return 1./(1. + exp(-x));
}

SEXP inv_logit_mat(SEXP x, SEXP min, SEXP max)
{
  int i, j, ind;
  const int m = nrows(x);
  const int n = ncols(x);
  double xval, pval;
  SEXP p;
  PROTECT(p = allocMatrix(REALSXP, m, n));
  
  
  for (j=0; j<n; j++)
  {
    for (i=0; i<m; i++)
    {
      ind = i + m*j;
      xval = REAL(x)[ind];
      pval = logistic(xval);
      
      if (ISNA(pval) && !ISNA(xval))
        pval = 1.;
      
      REAL(p)[ind] = pval * (REAL(max)[0] - REAL(min)[0]) + REAL(min)[0];
    }
  }
  
  UNPROTECT(1);
  return p;
}



/*this_loglike = sum(log(inv.logit.mat(q * theta))[q != 0])*/
SEXP compute_loglik(SEXP q, SEXP theta)
{
  double loglik = 0;
  int i, j, ind;
  const int m = nrows(q);
  const int n = ncols(q);
  double qval, tmp;
  SEXP ret;
  PROTECT(ret = allocVector(REALSXP, 1));
  
  
  for (j=0; j<n; j++)
  {
    for (i=0; i<m; i++)
    {
      ind = i + m*j;
      qval = REAL(q)[ind] * REAL(theta)[ind];
      tmp = logistic(qval);
      
      if (ISNA(tmp) && !ISNA(qval))
        tmp = 1.;
      
      loglik += log(tmp);
    }
  }
  
  REAL(ret)[0] = loglik;
  
  UNPROTECT(1);
  return ret;
}
















