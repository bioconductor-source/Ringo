#include <vector>
#include <algorithm>

using namespace std;

// the above stuff is pure C++, thus needs to be outside the 'extern "C"'

extern "C" {

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include <R_ext/Rdynload.h>
 
  //extern "C" {
  SEXP sliding_median(SEXP, SEXP, SEXP);
  //};

  /* HEADER: LOGISTICS MOSTLY IMPORTANT FOR WINDOWS DLL */

  /* this part just does not work so far, maybe needs to be in a .c file 
     or within extern "C"*/

  static R_CallMethodDef Ringo_calls[] = {
    {"sliding_median", (DL_FUNC) &sliding_median, 3},

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


  /* This function compute the median within a sliding window.
     Since a irregular spacing of data points is allowed, the number of 
     data points inside the window varies with the window position.
     Therefore the number of points inside the window is returned in addition
     to the median value.

     Written by Oleg Sklyar and Joern Toedling, January 2007.
  */

  SEXP
  sliding_median(SEXP ind, SEXP val, SEXP hlfsize) {
    int * x, nval, hs, is;
    double * y, * rval, median;
    int nprotect = 0;

    /* for the list of values of which we take the median at each position 
       in 'ind', we use a C++ 'vector' object, since we need to dynamically
       resize it */

    vector<double> values;

    SEXP res, dim;

    x = INTEGER(ind);
    y = REAL(val);

    hs = INTEGER(hlfsize)[0];
    nval = LENGTH(ind);

    PROTECT( res = allocVector( REALSXP, nval * 2) );
    nprotect++;

    rval = REAL(res);

    int i, j;

    for (i = 0; i < nval * 2; i++)
      rval[i] = R_NaN;

    int first = 0, last = 0;
    // we assume indexes, and therefore x, to be in accending order!
    // iterate over all positions, center the sliding window at position i:
    for (i = 0; i < nval; i++ ) {

      /* see if due to moving the window, values at left side are dropped
	 can be dropped window since they are more than half.width hs distant: */
      for (; first <= last && x[first] < x[i] - hs; ++first);

      /* see if due to moving window, values at right side should be appended
	 can be dropped window since they are more than half.width hs distant: */
      while ( x[last + 1] <= x[i] + hs && last < nval - 1) last++;
      
      is = last - first + 1; 
      // number of positions over which the median is taken

      if ( is == 0 ) continue; // if something went wrong
      values.resize(is); // 'resize' is a method of vector objects

      for (j = 0; j < is; j++)
	values[j] = y[first + j]; // get the values considered for median

      sort(values.begin(), values.end()); // sort values
      
      /* now median of sorted list of odd length is defined as the element in
	 in the middle element. If list has even length it's the mean of the
	 the two elements in the middle */
      j = is / 2;
      if ( j * 2 != is )
	rval[i] = values[j];
      else
	rval[i] = (values[j] + values[j-1]) / 2.0;
      rval[i + nval] = is;
    }
    // set the dimensions of the result: it's an array with number of
    //  positions rows and two columns
    PROTECT( dim = allocVector(INTSXP, 2));
    nprotect++;
    INTEGER(dim)[0] = nval;
    INTEGER(dim)[1] = 2;
    SET_DIM(res, dim);

    // do not need to protect objects from R's garbage collector any longer
    UNPROTECT(nprotect);

    return res;
  } 

}
