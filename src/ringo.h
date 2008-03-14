#include <R.h>
#include <Rdefines.h>

#ifdef __cplusplus
extern "C" {
#endif

SEXP sliding_quantile(SEXP, SEXP, SEXP, SEXP);
SEXP moving_mean_sd(SEXP, SEXP, SEXP);
SEXP overlap_xy(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

#ifdef __cplusplus
};
#endif
