#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _gicf_gicf_wrapper(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gicf_profileloglik(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_gicf_gicf_wrapper",  (DL_FUNC) &_gicf_gicf_wrapper,  10},
  {"_gicf_profileloglik", (DL_FUNC) &_gicf_profileloglik,  3},
  {NULL, NULL, 0}
};

void R_init_gicf(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
