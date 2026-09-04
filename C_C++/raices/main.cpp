#include <iostream>
#include <iomanip>
#include <cmath>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

double f(double x, void *params) {
  return 0;
}

int main() {
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  gsl_function F;
  F.function = &f;
  F.params = nullptr;
  double x_lo = 0.0;
  double x_hi = 1.0;
  T = gsl_root_fsolver_bisection;
  s = gsl_root_fsolver_alloc(T);
  gsl_root_fsolver_set(s, &F, x_lo, x_hi);
  std::cout << "iter\t" << "inf\t" << "sup\t" << "raíz\n";

  int status;
  int iter = 0;
  int max_iter = 100;
  double r;

  do {
    iter++;
    status = gsl_root_fsolver_iterate(s);
    r = gsl_root_fsolver_root(s);
    x_lo = gsl_root_fsolver_x_lower(s);
    x_hi = gsl_root_fsolver_x_upper(s);
    std::cout << iter << "\t" << x_lo << "\t" << x_hi << "\t" << r << "\n"; 
          status = gsl_root_test_interval( x_lo, x_hi, 0.0, 1e-8);
  } while(status == GSL_CONTINUE && iter < max_iter);

  std::cout << "\nRaiz encontrada = " << r << std::endl; gsl_root_fsolver_free(s);
  return 0;
}

