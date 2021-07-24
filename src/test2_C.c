#include <R.h>
#include <Rinternals.h>
#include <omp.h>

SEXP testeC(SEXP x_, SEXP y_) {
  double x = asReal(x_);
  double y = asReal(y_);

  double sum = x + y;

  return ScalarReal(sum);
}
