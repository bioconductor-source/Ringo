// initialization of the package

#include "ringo.h"

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/RConverters.h>
#include <R_ext/Rdynload.h>
//#include <R_ext/Utils.h>

/* LOGISTICS MOSTLY IMPORTANT FOR WINDOWS DLL */

static R_CallMethodDef Ringo_calls[] = {
  {"sliding_quantile", (DL_FUNC) &sliding_quantile, 4},
  {"moving_mean_sd", (DL_FUNC) &moving_mean_sd, 3},
  {"overlap_xy", (DL_FUNC) &overlap_xy, 6},
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
