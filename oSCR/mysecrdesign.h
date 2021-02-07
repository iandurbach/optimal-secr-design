/* constants */
  
  #define huge 1e10
  
  /* source to include */
  #include <math.h>
  #include <stdlib.h>
  #include <stdio.h>
  #include <R.h>       /* random numbers */
  #include <Rmath.h>   /* R math functions e.g. dbinom, dpois */
  #include <R_ext/Applic.h>
  #include <R_ext/Utils.h>
  
  void myLambdaC (
    double *par,       /* lambda0, sigma, z */
      int    *kk,        /* number of traps */
      int    *mm,        /* number of points on mask */
      double *traps,     /* x,y locations of traps (first x, then y) */
      double *mask,      /* x,y points on mask (first x, then y) */
      int    *fn,        /* detectfn code 0 = halfnormal */
      double *L,         /* return value vector of length mm */
      int    *resultcode /* 0 for successful completion */
  );
void mysumpkC (
  int    *type,      /* 0 multi 1 proximity */
    double *par,       /* lambda0, sigma, z */
    int    *kk,        /* number of traps */
    int    *mm,        /* number of points on mask */
    double *traps,     /* x,y locations of traps (first x, then y) */
    double *mask,      /* x,y points on mask (first x, then y) */
    int    *fn,        /* detectfn code 0 = halfnormal */
    double *L,         /* return value vector of length mm */
    double *L2,        /* return value vector of length mm */
    int    *resultcode /* 0 for successful completion */
);