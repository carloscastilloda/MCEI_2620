// Configuración de Librerias
#include <iostream>
#include <iomanip>
#include <cmath>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

// Función de trabajo: f(x) = x^3 - 5x + 1
double f(double x, void *params) {
  (void)params;
  return x * x * x - 5.0 * x + 1.0;
}
// Derivada de f(x) = x^3 - 5x + 1 es f'(x) = 3x^2 - 5
double df(double x, void *params) {
  (void)params;
  return 3.0 * x * x - 5.0;
}

// comienzo del código
// Funcion que busca una raíz dentro del intervalo [x_lo, x_hi] por metodo de biseccion
int buscar_raiz_biseccion(double x_lo, double x_hi) {
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  gsl_function F;

  F.function = &f;
  F.params   = nullptr;

  // Verificación previa del cambio de signo (Teorema de Bolzano)
  if (f(x_lo, nullptr) * f(x_hi, nullptr) > 0.0) {
    std::cout << "  [!] No hay cambio de signo en [" << x_lo << ", " << x_hi << "]\n";
    return -1;
  }

  // uso de gsl para resolver raices por biseccion
  T = gsl_root_fsolver_bisection;
  s = gsl_root_fsolver_alloc(T);
  gsl_root_fsolver_set(s, &F, x_lo, x_hi);

  std::cout << "  Solver: " << gsl_root_fsolver_name(s) << "   intervalo inicial ["
            << x_lo << ", " << x_hi << "]\n";
  std::cout << "  " << std::setw(5) << "iter" << std::setw(16) << "inf"
            << std::setw(16) << "sup" << std::setw(16) << "raiz"
            << std::setw(14) << "amplitud" << "\n";
  std::cout << "  " << std::string(67, '-') << "\n";

  // asignacion de variables status, iter, max_iter y r
  int status;
  int iter = 0;
  int max_iter = 100;
  double r = 0.0;

  std::cout << std::fixed << std::setprecision(10);

  do {
    iter++;
    status = gsl_root_fsolver_iterate(s);
    r      = gsl_root_fsolver_root(s);
    x_lo   = gsl_root_fsolver_x_lower(s);
    x_hi   = gsl_root_fsolver_x_upper(s);

    std::cout << "  " << std::setw(5) << iter << std::setw(16) << x_lo
              << std::setw(16) << x_hi << std::setw(16) << r
              << std::setw(14) << (x_hi - x_lo) << "\n";

    // criterio de parada: amplitud del intervalo < 1e-8
    status = gsl_root_test_interval(x_lo, x_hi, 0.0, 1e-8);
  } while (status == GSL_CONTINUE && iter < max_iter);

  std::cout << "  ->  Raiz = " << r
            << "   |  f(raiz) = " << std::scientific << f(r, nullptr)
            << std::fixed << "   |  iteraciones = " << iter
            << (status == GSL_SUCCESS ? "   [CONVERGE]\n" : "   [NO CONVERGE]\n");

  gsl_root_fsolver_free(s);
  return status;
}

// Funcion buscar raiz por metodo de falsa posición//
int buscar_raiz_falsa_posicion(double x_lo, double x_hi) {
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  gsl_function F;

  F.function = &f;
  F.params   = nullptr;

  if (f(x_lo, nullptr) * f(x_hi, nullptr) > 0.0) {
    std::cout << "  [!] No hay cambio de signo en [" << x_lo << ", " << x_hi << "]\n";
    return -1;
  }

  // uso de gsl para resolver raices por falsa posicion //
  T = gsl_root_fsolver_falsepos;
  s = gsl_root_fsolver_alloc(T);
  gsl_root_fsolver_set(s, &F, x_lo, x_hi);

  std::cout << "  Solver: " << gsl_root_fsolver_name(s) << "   intervalo inicial ["
            << x_lo << ", " << x_hi << "]\n";
  std::cout << "  " << std::setw(5) << "iter" << std::setw(16) << "inf"
            << std::setw(16) << "sup" << std::setw(16) << "raiz"
            << std::setw(14) << "amplitud" << "\n";
  std::cout << "  " << std::string(67, '-') << "\n";

  int status;
  int iter = 0;
  int max_iter = 100;
  double r = 0.0;

  std::cout << std::fixed << std::setprecision(10);
  do {
    iter++;
    status = gsl_root_fsolver_iterate(s);
    r      = gsl_root_fsolver_root(s);
    x_lo   = gsl_root_fsolver_x_lower(s);
    x_hi   = gsl_root_fsolver_x_upper(s);

    std::cout << "  " << std::setw(5) << iter << std::setw(16) << x_lo
              << std::setw(16) << x_hi << std::setw(16) << r
              << std::setw(14) << (x_hi - x_lo) << "\n";

    status = gsl_root_test_interval(x_lo, x_hi, 0.0, 1e-8);
  } while (status == GSL_CONTINUE && iter < max_iter);

  std::cout << "  ->  Raiz = " << r
            << "   |  f(raiz) = " << std::scientific << f(r, nullptr)
            << std::fixed << "   |  iteraciones = " << iter
            << (status == GSL_SUCCESS ? "   [CONVERGE]\n" : "   [NO CONVERGE]\n");

  gsl_root_fsolver_free(s);
  return status;
}

// Funcion buscar raiz por metodo de brent//
int buscar_raiz_brent(double x_lo, double x_hi) {
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  gsl_function F;

  F.function = &f;
  F.params   = nullptr;

  if (f(x_lo, nullptr) * f(x_hi, nullptr) > 0.0) {
    std::cout << "  [!] No hay cambio de signo en [" << x_lo << ", " << x_hi << "]\n";
    return -1;
  }

  // uso de gsl para resolver raices por el metodo de Brent
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc(T);
  gsl_root_fsolver_set(s, &F, x_lo, x_hi);

  std::cout << "  Solver: " << gsl_root_fsolver_name(s) << "   intervalo inicial ["
            << x_lo << ", " << x_hi << "]\n";
  std::cout << "  " << std::setw(5) << "iter" << std::setw(16) << "inf"
            << std::setw(16) << "sup" << std::setw(16) << "raiz"
            << std::setw(14) << "amplitud" << "\n";
  std::cout << "  " << std::string(67, '-') << "\n";

  int status;
  int iter = 0;
  int max_iter = 100;
  double r = 0.0;

  std::cout << std::fixed << std::setprecision(10);

  do {
    iter++;
    status = gsl_root_fsolver_iterate(s);
    r      = gsl_root_fsolver_root(s);
    x_lo   = gsl_root_fsolver_x_lower(s);
    x_hi   = gsl_root_fsolver_x_upper(s);

    std::cout << "  " << std::setw(5) << iter << std::setw(16) << x_lo
              << std::setw(16) << x_hi << std::setw(16) << r
              << std::setw(14) << (x_hi - x_lo) << "\n";

    status = gsl_root_test_interval(x_lo, x_hi, 0.0, 1e-8);
  } while (status == GSL_CONTINUE && iter < max_iter);

  std::cout << "  ->  Raiz = " << r
            << "   |  f(raiz) = " << std::scientific << f(r, nullptr)
            << std::fixed << "   |  iteraciones = " << iter
            << (status == GSL_SUCCESS ? "   [CONVERGE]\n" : "   [NO CONVERGE]\n");

  gsl_root_fsolver_free(s);
  return status;
}

// Cálculo de raices por newton 
// GSL pide una rutina que entregue f y f' a la vez
void fdf(double x, void *params, double *y, double *dy) {
  (void)params;
  *y  = x * x * x - 5.0 * x + 1.0;
  *dy = 3.0 * x * x - 5.0;
}

int buscar_raiz_newton(double x) {
  const gsl_root_fdfsolver_type *T;
  gsl_root_fdfsolver *s;
  gsl_function_fdf FDF;

  FDF.f      = &f;
  FDF.df     = &df;
  FDF.fdf    = &fdf;
  FDF.params = nullptr;

  // Aviso: si f'(x0) = 0 el metodo falla (division por cero)
  if (std::fabs(df(x, nullptr)) < 1e-12) {
    std::cout << "  [!] f'(x0) ~ 0, el metodo de Newton no puede arrancar en x0 = "
              << x << "\n";
    return -1;
  }

  // uso de gsl para resolver raices por Newton
  T = gsl_root_fdfsolver_newton;
  s = gsl_root_fdfsolver_alloc(T);
  gsl_root_fdfsolver_set(s, &FDF, x);

  std::cout << "  Solver: " << gsl_root_fdfsolver_name(s)
            << "   valor inicial x0 = " << x << "\n";
  std::cout << "  " << std::setw(5) << "iter" << std::setw(18) << "x_anterior"
            << std::setw(18) << "raiz" << std::setw(16) << "|dx|" << "\n";
  std::cout << "  " << std::string(57, '-') << "\n";

  // asignacion de variables status, iter, max_iter y x0
  int status;
  int iter = 0;
  int max_iter = 100;
  double x0;

  std::cout << std::fixed << std::setprecision(10);

  // ejecutar el calculo de la raiz por Newton
  do {
    iter++;
    x0     = x;                                  // se guarda la aproximacion previa
    status = gsl_root_fdfsolver_iterate(s);
    if (status != GSL_SUCCESS) break;            // p.ej. derivada nula
    x      = gsl_root_fdfsolver_root(s);

    std::cout << "  " << std::setw(5) << iter << std::setw(18) << x0
              << std::setw(18) << x << std::setw(16) << std::fabs(x - x0) << "\n";

    // criterio de parada: |x_n - x_{n-1}| < 1e-8
    status = gsl_root_test_delta(x, x0, 0.0, 1e-8);
  } while (status == GSL_CONTINUE && iter < max_iter);

  std::cout << "  ->  Raiz = " << x
            << "   |  f(raiz) = " << std::scientific << f(x, nullptr)
            << std::fixed << "   |  iteraciones = " << iter
            << (status == GSL_SUCCESS ? "   [CONVERGE]\n" : "   [NO CONVERGE]\n");

  gsl_root_fdfsolver_free(s);
  return status;
}

//Metodo Secante
int buscar_raiz_secante(double x) {
  const gsl_root_fdfsolver_type *T;
  gsl_root_fdfsolver *s;
  gsl_function_fdf FDF;

  FDF.f      = &f;
  FDF.df     = &df;
  FDF.fdf    = &fdf;
  FDF.params = nullptr;

  // uso de gsl para resolver raices por el metodo de la secante
  T = gsl_root_fdfsolver_secant;
  s = gsl_root_fdfsolver_alloc(T);
  gsl_root_fdfsolver_set(s, &FDF, x);

  std::cout << "  Solver: " << gsl_root_fdfsolver_name(s)
            << "   valor inicial x0 = " << x << "\n";
  std::cout << "  " << std::setw(5) << "iter" << std::setw(18) << "x_anterior"
            << std::setw(18) << "raiz" << std::setw(16) << "|dx|" << "\n";
  std::cout << "  " << std::string(57, '-') << "\n";

  int status;
  int iter = 0;
  int max_iter = 100;
  double x0;

  std::cout << std::fixed << std::setprecision(10);

  do {
    iter++;
    x0     = x;
    status = gsl_root_fdfsolver_iterate(s);
    if (status != GSL_SUCCESS) break;
    x      = gsl_root_fdfsolver_root(s);

    std::cout << "  " << std::setw(5) << iter << std::setw(18) << x0
              << std::setw(18) << x << std::setw(16) << std::fabs(x - x0) << "\n";

    status = gsl_root_test_delta(x, x0, 0.0, 1e-8);
  } while (status == GSL_CONTINUE && iter < max_iter);

  std::cout << "  ->  Raiz = " << x
            << "   |  f(raiz) = " << std::scientific << f(x, nullptr)
            << std::fixed << "   |  iteraciones = " << iter
            << (status == GSL_SUCCESS ? "   [CONVERGE]\n" : "   [NO CONVERGE]\n");

  gsl_root_fdfsolver_free(s);
  return status;
}

//Raices por Método Steffenson
int buscar_raiz_steffenson(double x) {
  const gsl_root_fdfsolver_type *T;
  gsl_root_fdfsolver *s;
  gsl_function_fdf FDF;

  FDF.f      = &f;
  FDF.df     = &df;
  FDF.fdf    = &fdf;
  FDF.params = nullptr;

  if (std::fabs(df(x, nullptr)) < 1e-12) {
    std::cout << "  [!] f'(x0) ~ 0, Steffenson no puede arrancar en x0 = " << x << "\n";
    return -1;
  }

  // uso de gsl para resolver raices por el metodo de Steffenson
  T = gsl_root_fdfsolver_steffenson;
  s = gsl_root_fdfsolver_alloc(T);
  gsl_root_fdfsolver_set(s, &FDF, x);

  std::cout << "  Solver: " << gsl_root_fdfsolver_name(s)
            << "   valor inicial x0 = " << x << "\n";
  std::cout << "  " << std::setw(5) << "iter" << std::setw(18) << "x_anterior"
            << std::setw(18) << "raiz" << std::setw(16) << "|dx|" << "\n";
  std::cout << "  " << std::string(57, '-') << "\n";

  int status;
  int iter = 0;
  int max_iter = 100;
  double x0;

  std::cout << std::fixed << std::setprecision(10);

  do {
    iter++;
    x0     = x;
    status = gsl_root_fdfsolver_iterate(s);
    if (status != GSL_SUCCESS) {
      std::cout << "  [!] " << gsl_strerror(status) << "\n";
      break;
    }
    x = gsl_root_fdfsolver_root(s);

    std::cout << "  " << std::setw(5) << iter << std::setw(18) << x0
              << std::setw(18) << x << std::setw(16) << std::fabs(x - x0) << "\n";

    status = gsl_root_test_delta(x, x0, 0.0, 1e-8);
  } while (status == GSL_CONTINUE && iter < max_iter);

  std::cout << "  ->  Raiz = " << x
            << "   |  f(raiz) = " << std::scientific << f(x, nullptr)
            << std::fixed << "   |  iteraciones = " << iter
            << (status == GSL_SUCCESS ? "   [CONVERGE]\n" : "   [NO CONVERGE]\n");

  gsl_root_fdfsolver_free(s);
  return status;
}

// Función Principal
int main() {
  const double intervalos[3][2] = {
      {-3.0, -2.0},
      { 0.0,  1.0},
      { 2.0,  3.0}
  };

  const double x0[3] = {-2.5, 0.0, 2.5};

  std::cout << " Función -  f(x) = x^3 - 5x + 1\n";
  std::cout << " Derivada de la Función -  f(x) = 3x^2 - 5\n";

  for (int k = 0; k < 3; ++k) {
    std::cout << "\n>> Raiz Biseccion #" << (k + 1) << "\n";
    buscar_raiz_biseccion(intervalos[k][0], intervalos[k][1]);
    std::cout << "\n>> Raiz Falsa Posición#" << (k + 1) << "\n";
    buscar_raiz_falsa_posicion(intervalos[k][0], intervalos[k][1]);
    std::cout << "\n>> Raiz Método Brent#" << (k + 1) << "\n";
    buscar_raiz_brent(intervalos[k][0], intervalos[k][1]);
    std::cout << "\n>> Raiz Método Newton #" << (k + 1) << "\n";
    buscar_raiz_newton(x0[k]);
    std::cout << "\n>> Raiz Método Secante#" << (k + 1) << "\n";
    buscar_raiz_secante(x0[k]);
    std::cout << "\n>> Raiz Método Steffenson#" << (k + 1) << "\n";
    buscar_raiz_steffenson(x0[k]);
  }

  return 0;
}