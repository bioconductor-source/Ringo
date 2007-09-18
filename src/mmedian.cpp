#include <vector>
#include <algorithm>
#include <list>

using namespace std;

// the above stuff is pure C++, thus needs to be outside the 'extern "C"'

extern "C" {

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Error.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Utils.h> 
 
  //extern "C" {
  SEXP sliding_median(SEXP, SEXP, SEXP);
  SEXP sliding_quantile(SEXP, SEXP, SEXP, SEXP);
  //};

  /* HEADER: LOGISTICS MOSTLY IMPORTANT FOR WINDOWS DLL */

   static R_CallMethodDef Ringo_calls[] = {
     //{"sliding_median", (DL_FUNC) &sliding_median, 3},
    {"sliding_quantile", (DL_FUNC) &sliding_quantile, 4},

    // necessary last entry of R_CallMethodDef:
    {NULL, NULL, 0}
  };

  // register functions 
  void R_init_Ringo(DllInfo *dll)
  {
    R_registerRoutines(dll, NULL, Ringo_calls, NULL, NULL);
  }

  void R_unload_Ringo(DllInfo *dll) 
  {
    // at the moment nothing to do here 
  }

  /* MAIN PART */


 /* This next function replaced the sliding.median function
     above. It uses a different data structure, a linked list, and in 
     most use cases it is faster than the sliding.median function.
     It also allows the computation of any quantile in the window, not
     only the median. Again, irregular spacing of data points.
     The number of points inside the window is returned in addition
     to the quantile value.

     Written by Oleg Sklyar and Joern Toedling, March 2007.
  */

  SEXP
  sliding_quantile(SEXP ind, SEXP val, SEXP hlfsize, SEXP probability) {
    /* ind and x must be passed sorted! ensured in R */
    int i, j, nprotect = 0;

    SEXP res;

    int * x = INTEGER(ind);
    double * y = REAL(val);

    int hs = INTEGER(hlfsize)[0];
    int nval = LENGTH(ind);
    double prob = REAL(probability)[0];

    /* result, initialized with NA */
    PROTECT( res = allocVector( REALSXP, nval * 2) );
    nprotect++;
    double * rval = REAL(res);
    for (i = 0; i < nval * 2; i++) rval[i] = R_NaN;

    /* we keep two lists, yy always sorted list of y-values currently in
     * the frame and xx - a list of corresponding indexes. Every time the
     * frame moves we calculate the 'first' and 'last' indexes of points
     * in the original vectors ind and val that currently belong to this
     * list. Then we drop all those from the lists, for which indexes are
     * below first and we add (in correct positions to preserve sorting)
     * the values between previous last (excl) and current last (incl) */
    list<int> xx;
    list<int>::iterator xit, tmpxit;
    list<double> yy;
    list<double>::iterator yit, tmpyit;

    int first = 0, last = -1, prevlast, pos;

    for (i = 0; i < nval; i++ ) { // i is position of the frame middle
      /* calculate first index */
      for (; first <= last && x[first] < x[i] - hs; ++first);
      /* drop values from the list that have indexes below first */
      if ( xx.size() > 0 ) {
        yit = yy.begin();
        xit = xx.begin();
        while ( xit != xx.end() ) {
          if ( *xit < first ) {
            tmpxit = xit; xit++;
            tmpyit = yit; yit++;
            xx.erase(tmpxit);
            yy.erase(tmpyit);
            continue; // moves to next loop iteration, like "next" in R
          }
          xit++;
          yit++;
        }
      }

      /* backup current last */
      prevlast = last;
      /* calculate new last */
      while ( x[last + 1] <= x[i] + hs && last < nval - 1) last++;
      /* add values between prevlast+1 and last to the list, preserving sort */

      xit = xx.begin();
      yit = yy.begin();
      for ( j = prevlast + 1; j <= last; j++ ) {
        if ( xx.size() == 0 ) {
          /* add first element to an empty list */
          xx.push_back( j );
          yy.push_back( y[j] );
          xit = xx.begin();
          yit = yy.begin();
          continue; // moves to next loop iteration, as "next" in R
        }
        // push it off the end
        if ( yit == yy.end() ) {
          yit--;
          xit--;
        }
        // decrease iterator until *yit becomes <= y[j]
        while ( yit != yy.begin() && *yit >= y[j] ) {
          yit--;
          xit--;
        }
        // then move one step upwards or until we get over if we started low
        while ( yit != yy.end() && *yit < y[j] ) {
          yit++;
          xit++;
        }
        xx.insert(xit, j );
        yy.insert(yit, y[j]);
      }

      if ( yy.size() == 0 ) {
        UNPROTECT( nprotect ); // let the garbage collector delete res, which is not returned
        error("zero frame size in the sliding.quantile C-routine");
      }

      /* now the q (0<=q<=1) quantile of sorted list yy of length n is defined 
       * as follows:  first compute  k = (q * n)
       * then the quantile is
       * yy[k] if k is integer
       * yy[floor(k)] + q*(yy[ceiling(k)]-yy[floor(k)]),
       * the value of the element floor(k) plus the q proportion of the 
       * difference between it and the next higher element */

      pos = (int)((yy.size() - 1) * prob);
      yit = yy.begin();
      for ( j = 0; j < pos; j++ ) yit++;
      rval[i] = *yit;
      if ( (int)( pos / prob ) + 1 != (int) yy.size() ) {
        yit++;
        rval[i] = *yit * prob + rval[i] * (1.0 - prob);
      }
      rval[i + nval] = yy.size();
      R_CheckUserInterrupt();
    }/*for i where positions[i] is the middle position of the sliding window*/
   
    /* set the dimensions of the result: it's an array with number of
     *  positions rows and two columns */
    SEXP dim;
    PROTECT( dim = allocVector(INTSXP, 2));
    nprotect++;
    INTEGER(dim)[0] = nval;
    INTEGER(dim)[1] = 2;
    SET_DIM(res, dim);

    UNPROTECT(nprotect);
    return res;
  } // sliding_quantile

}// extern C
