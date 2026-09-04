// Configuración de Librerias
#include <iostream>
#include <iomanip>
#include <cmath>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

// Identificador de la funcion de trabajo (se pasa por *params)
enum Funcion { CUBICA = 0, EXPONENCIAL = 1 };
 
// Funcion de trabajo: se elige segun el valor apuntado por params
double f(double x, void *params) {
  const Funcion id = *static_cast<Funcion *>(params);
  switch (id) {
    case CUBICA:      return x * x * x - 5.0 * x + 1.0;
    case EXPONENCIAL: return std::exp(x) - x;
  }
  return 0.0;
}

// Derivada de la funcion de trabajo
double df(double x, void *params) {
  const Funcion id = *static_cast<Funcion *>(params);
  switch (id) {
    case CUBICA:      return 3.0 * x * x - 5.0;
    case EXPONENCIAL: return std::exp(x) - 1.0;
  }
  return 0.0;
}

// Son exactamente 4 parametros: uno para f y uno para f', no dos de cada.
void fdf(double x, void *params, double *y, double *dy) {
  *y  = f(x, params);
  *dy = df(x, params);
}

// Descripcion de cada problema: funcion, intervalos y semillas propios
struct Problema {
  Funcion     id;
  const char *nombre;
  const char *derivada;
  int         n_casos;
  double      intervalos[3][2];   // para los metodos cerrados
  double      semillas[3];        // para los metodos abiertos
  const char *nota;
};

// METODOS CERRADOS  (gsl_root_fsolver: bisección, falsa posición, Brent)
int buscar_raiz_cerrado(const gsl_root_fsolver_type *T, Funcion id,
                        double x_lo, double x_hi) {
  gsl_function F;
  F.function = &f;
  F.params   = &id;                 // <-- aqui se selecciona la funcion
 
  // Verificación previa del cambio de signo (Teorema de Bolzano)
  if (f(x_lo, &id) * f(x_hi, &id) > 0.0) {
    std::cout << "  [!] No hay cambio de signo en [" << x_lo << ", " << x_hi
              << "]  ->  el metodo cerrado no es aplicable\n";
    return -1;
  }
 
  gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
  gsl_root_fsolver_set(s, &F, x_lo, x_hi);
 
  std::cout << "  Solver: " << gsl_root_fsolver_name(s) << "   intervalo inicial ["
            << x_lo << ", " << x_hi << "]\n";
  std::cout << "  " << std::setw(5) << "iter" << std::setw(16) << "inf"
            << std::setw(16) << "sup" << std::setw(16) << "raiz"
            << std::setw(14) << "amplitud" << "\n";
  std::cout << "  " << std::string(67, '-') << "\n";
 
  // asignacion de variables status, iter, max_iter y r
  int status = GSL_CONTINUE;
  int iter = 0;
  int max_iter = 100;
  double r = 0.0;
 
  std::cout << std::fixed << std::setprecision(10);
 
  do {
    iter++;
    if (gsl_root_fsolver_iterate(s) != GSL_SUCCESS) break;
    r    = gsl_root_fsolver_root(s);
    x_lo = gsl_root_fsolver_x_lower(s);
    x_hi = gsl_root_fsolver_x_upper(s);
 
    std::cout << "  " << std::setw(5) << iter << std::setw(16) << x_lo
              << std::setw(16) << x_hi << std::setw(16) << r
              << std::setw(14) << (x_hi - x_lo) << "\n";
 
    // criterio de parada: amplitud del intervalo < 1e-8
    status = gsl_root_test_interval(x_lo, x_hi, 0.0, 1e-8);
  } while (status == GSL_CONTINUE && iter < max_iter);
 
  std::cout << "  ->  Raiz = " << r
            << "   |  f(raiz) = " << std::scientific << f(r, &id)
            << std::fixed << "   |  iteraciones = " << iter
            << (status == GSL_SUCCESS ? "   [CONVERGE]\n" : "   [NO CONVERGE]\n");
 
  gsl_root_fsolver_free(s);
  return status;
}

// METODOS ABIERTOS  (gsl_root_fdfsolver: Newton, secante, Steffenson)
int buscar_raiz_abierto(const gsl_root_fdfsolver_type *T, Funcion id, double x) {
  gsl_function_fdf FDF;
  FDF.f      = &f;
  FDF.df     = &df;
  FDF.fdf    = &fdf;
  FDF.params = &id;                 // Selección de la funcion
 
  // si f'(x0) = 0 el metodo no funciona indeterminado / division por cero
  if (std::fabs(df(x, &id)) < 1e-12) {
    std::cout << "  [!] f'(x0) ~ 0, el metodo no puede arrancar en x0 = " << x << "\n";
    return -1;
  }
 
  gsl_root_fdfsolver *s = gsl_root_fdfsolver_alloc(T);
  gsl_root_fdfsolver_set(s, &FDF, x);
 
  std::cout << "  Solver: " << gsl_root_fdfsolver_name(s)
            << "   valor inicial x0 = " << x << "\n";
  std::cout << "  " << std::setw(5) << "iter" << std::setw(18) << "x_anterior"
            << std::setw(18) << "raiz" << std::setw(16) << "|dx|" << "\n";
  std::cout << "  " << std::string(57, '-') << "\n";
 
  // asignacion de variables status, iter, max_iter y x0
  int status = GSL_CONTINUE;
  int iter = 0;
  int max_iter = 100;
  double x0;
 
  std::cout << std::fixed << std::setprecision(10);
 
  do {
    iter++;
    x0 = x;                                      // se guarda la aproximacion previa
    int it_status = gsl_root_fdfsolver_iterate(s);
    if (it_status != GSL_SUCCESS) {             
      std::cout << "  [!] " << gsl_strerror(it_status) << "\n";
      status = it_status;
      break;
    }
    x = gsl_root_fdfsolver_root(s);
 
    std::cout << "  " << std::setw(5) << iter << std::setw(18) << x0
              << std::setw(18) << x << std::setw(16) << std::fabs(x - x0) << "\n";
 
    // criterio de parada: |x_n - x_{n-1}| < 1e-8
    status = gsl_root_test_delta(x, x0, 0.0, 1e-8);
  } while (status == GSL_CONTINUE && iter < max_iter);
 
  const double residuo = std::fabs(f(x, &id));
  const bool   es_raiz = (residuo < 1e-6);
 
  std::cout << "  ->  Raiz = " << x
            << "   |  f(raiz) = " << std::scientific << f(x, &id)
            << std::fixed << "   |  iteraciones = " << iter;
  if (status == GSL_SUCCESS && es_raiz)
    std::cout << "   [CONVERGE]\n";
  else if (status == GSL_SUCCESS)
    std::cout << "   [FALSO POSITIVO: |dx| pequeno pero f(x) != 0]\n";
  else
    std::cout << "   [NO CONVERGE]\n";
 
  gsl_root_fdfsolver_free(s);
  return (status == GSL_SUCCESS && es_raiz) ? GSL_SUCCESS : GSL_CONTINUE;
}

// Ejecución de los 6 metodos sobre la ecuación del problema
void resolver(const Problema &p) {
  std::cout << " FUNCION:  " << p.nombre << "\n"
            << " DERIVADA: " << p.derivada << "\n";
  if (p.nota) std::cout << " " << p.nota << "\n";
 
  for (int k = 0; k < p.n_casos; ++k) {
    std::cout << "Caso " << (k + 1) << " de " << p.n_casos
              << "   intervalo [" << p.intervalos[k][0] << ", " << p.intervalos[k][1]
              << "]   semilla x0 = " << p.semillas[k] << "/n";
 
    std::cout << "\n>> 1. Biseccion\n";
    buscar_raiz_cerrado(gsl_root_fsolver_bisection, p.id,
                        p.intervalos[k][0], p.intervalos[k][1]);
 
    std::cout << "\n>> 2. Falsa posicion\n";
    buscar_raiz_cerrado(gsl_root_fsolver_falsepos, p.id,
                        p.intervalos[k][0], p.intervalos[k][1]);
 
    std::cout << "\n>> 3. Brent\n";
    buscar_raiz_cerrado(gsl_root_fsolver_brent, p.id,
                        p.intervalos[k][0], p.intervalos[k][1]);
 
    std::cout << "\n>> 4. Newton\n";
    buscar_raiz_abierto(gsl_root_fdfsolver_newton, p.id, p.semillas[k]);
 
    std::cout << "\n>> 5. Secante\n";
    buscar_raiz_abierto(gsl_root_fdfsolver_secant, p.id, p.semillas[k]);
 
    std::cout << "\n>> 6. Steffenson\n";
    buscar_raiz_abierto(gsl_root_fdfsolver_steffenson, p.id, p.semillas[k]);
  }
}

// Función Principal
int main() {
  gsl_set_error_handler_off();   // sin esto GSL llama a abort() ante un error
 
  const Problema cubica = {
      CUBICA,
      "f(x) = x^3 - 5x + 1",
      "f'(x) = 3x^2 - 5",
      3,
      {{-3.0, -2.0}, {0.0, 1.0}, {2.0, 3.0}},
      {-2.5, 0.0, 2.5},
      "Tres raices reales. Cuidado: f'(x) = 0 en x = +-sqrt(5/3) ~ +-1.2910."
  };
 
  const Problema exponencial = {
      EXPONENCIAL,
      "f(x) = e^x - x",
      "f'(x) = e^x - 1",
      2,
      {{-2.0, 0.0}, {0.0, 2.0}, {0.0, 0.0}},
      {-1.0, 1.0, 0.0},
      "AVISO: esta funcion NO tiene raices reales. Su minimo global esta en\n"
      " x = 0 (donde f'(x) = e^x - 1 = 0) y vale f(0) = 1 > 0, es decir e^x > x\n"
      " para todo x real. Cuando NO se cumple la hipotesis de Bolzano:\n"
      " los metodos cerrados rechazan el intervalo y los abiertos no convergen.\n"
  };
 
  resolver(cubica);
  resolver(exponencial);
 
  std::cout << " Para usar otra funcion exponencial que SI tenga raiz, cambie el caso\n"
            << " EXPONENCIAL en f() y df(). Opciones habituales:\n"
            << "   e^-x - x     f' = -e^-x - 1   raiz  0.5671432904        [0, 1]\n"
            << "   e^x - x - 2  f' = e^x - 1     raices -1.8414, 1.1462    [-2,-1] [1,2]\n"
            << "   e^x - 3x     f' = e^x - 3     raices  0.6191, 1.5121    [0, 1] [1, 2]\n";
 
  return 0;
}