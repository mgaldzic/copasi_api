/* Begin CVS Header
   $Source: /fs/turing/cvs/copasi_dev/copasi/odepack++/dfnorm.cpp,v $
   $Revision: 1.3 $
   $Name: Build-33 $
   $Author: shoops $
   $Date: 2006/06/20 13:19:11 $
   End CVS Header */

// Copyright � 2006 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.
//
// This C++ code is based on an f2c conversion of the Fortran
// library ODEPACK available at: http://www.netlib.org/odepack/

#include <math.h>

#include <algorithm>

#include "copasi.h"

#include "dfnorm.h"

/* DECK DFNORM */
double dfnorm_(C_INT *n, double *a, double *w)
{
  /* System generated locals */
  C_INT a_dim1, a_offset, i__1, i__2;
  double ret_val, d__1, d__2;

  /* Local variables */
  C_INT i__, j;
  double an, sum;

  /* ----------------------------------------------------------------------- */
  /* This function computes the norm of a full N by N matrix, */
  /* stored in the array A, that is consistent with the weighted max-norm */
  /* on vectors, with weights stored in the array W: */
  /*   DFNORM = MAX(i=1,...,N) (W(i) * Sum(j=1,...,N) ABS(a(i,j))/W(j)) */
  /* ----------------------------------------------------------------------- */
  /* Parameter adjustments */
  --w;
  a_dim1 = *n;
  a_offset = 1 + a_dim1;
  a -= a_offset;

  /* Function Body */
  an = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      sum = 0.;
      i__2 = *n;
      for (j = 1; j <= i__2; ++j)
        {
          /* L10: */
          sum += (d__1 = a[i__ + j * a_dim1], fabs(d__1)) / w[j];
        }
      /* Computing MAX */
      d__1 = an, d__2 = sum * w[i__];
      an = std::max(d__1, d__2);
      /* L20: */
    }
  ret_val = an;
  return ret_val;
  /* ----------------------- End of Function DFNORM ------------------------ */
} /* dfnorm_ */
