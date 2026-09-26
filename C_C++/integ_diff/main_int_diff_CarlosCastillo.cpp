/* =====================================================================
 * Taller integración y diferenciación numérica - C + GSL
 * MCEI_M - Escuela Colombiana de Ingeniería
 * Elaborador Por: Carlos Castillo
 * ===================================================================== */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

#define NMAX 1000

/*  temporizador de alta resolución  */
static double ahora(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

/*  Ejercicio 1: función a integrar  */
typedef struct { double a, amp, w; } Params;

static double f(double x, void *p) {
    Params *q = (Params *)p;
    return exp(-q->a * x) * (1.0 + q->amp * sin(q->w * x));
}

static double referencia(void) {
    const double a = 0.4, b = 3.0, L = 8.0;
    return (1.0 - exp(-a * L)) / a
         + 0.5 * (b - exp(-a * L) * (a * sin(b * L) + b * cos(b * L))) / (a * a + b * b);
}

/* Trapecio compuesto explícito: el programador controla TODO */
static double trapecio_func(double (*g)(double, void *), void *p,
                            double a, double b, int n) {
    double h = (b - a) / n, s = 0.5 * (g(a, p) + g(b, p));
    for (int i = 1; i < n; i++) s += g(a + i * h, p);
    return h * s;
}

/*  Ejercicio 2: utilidades para datos  */
static int leer_csv(const char *ruta, double *x, double *y, int nmax) {
    FILE *fp = fopen(ruta, "r");
    if (!fp) { perror(ruta); return -1; }
    char linea[256];
    if (!fgets(linea, sizeof linea, fp)) { fclose(fp); return -1; } /* encabezado */
    int n = 0;
    while (n < nmax && fgets(linea, sizeof linea, fp)) {
        if (sscanf(linea, "%lf,%lf", &x[n], &y[n]) == 2) n++;
        else { fprintf(stderr, "Línea %d inválida o con faltantes\n", n + 2); }
    }
    fclose(fp);
    return n;
}

static double trapecio_datos(const double *x, const double *y, int n) {
    double s = 0.0;
    for (int i = 0; i < n - 1; i++) s += 0.5 * (x[i + 1] - x[i]) * (y[i] + y[i + 1]);
    return s;
}

static double simpson13(const double *y, int npts, double h) { /* npts impar */
    double s = y[0] + y[npts - 1];
    for (int i = 1; i < npts - 1; i++) s += (i % 2 ? 4.0 : 2.0) * y[i];
    return h / 3.0 * s;
}

static double simpson38(const double *y, double h) {             /* 4 puntos */
    return 3.0 * h / 8.0 * (y[0] + 3 * y[1] + 3 * y[2] + y[3]);
}

int main(void) {
    gsl_set_error_handler_off();   /* manejamos los códigos de estado nosotros */

    /*  EJERCICIO 1  */
    Params par = {0.4, 0.5, 3.0};
    gsl_function F = { .function = &f, .params = &par };
    const double a = 0.0, b = 8.0, Iref = referencia();
    printf("=== EJERCICIO 1 ===\nReferencia analítica: %.15f\n\n", Iref);

    /* --- Trapecio explícito --- */
    int N[] = {10, 20, 50, 100, 500, 1000};
    printf("%6s %20s %12s %12s\n", "n", "Trapecio", "Error abs", "Tiempo [s]");
    for (int k = 0; k < 6; k++) {
        int rep = 20000; double I = 0, t0 = ahora();
        for (int r = 0; r < rep; r++) I = trapecio_func(f, &par, a, b, N[k]);
        double t = (ahora() - t0) / rep;
        printf("%6d %20.15f %12.3e %12.3e\n", N[k], I, fabs(I - Iref), t);
    }

    /* --- Cuadratura adaptativa QAG (Gauss–Kronrod) --- */
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(NMAX);
    double tols[] = {1e-4, 1e-6, 1e-10, 1e-13};
    printf("\nQAG (GSL_INTEG_GAUSS21), epsabs = 0\n");
    printf("%8s %20s %12s %12s %8s %6s %12s\n",
           "epsrel", "Resultado", "Err. estim.", "Err. real", "status", "subint", "Tiempo [s]");
    for (int k = 0; k < 4; k++) {
        double res = 0, err = 0; int st = 0; int rep = 20000;
        double t0 = ahora();
        for (int r = 0; r < rep; r++)
            st = gsl_integration_qag(&F, a, b, 0.0, tols[k], NMAX,
                                     GSL_INTEG_GAUSS21, w, &res, &err);
        double t = (ahora() - t0) / rep;
        printf("%8.0e %20.15f %12.3e %12.3e %8s %6zu %12.3e\n",
               tols[k], res, err, fabs(res - Iref), gsl_strerror(st), w->size, t);
    }

    /* Efecto de la regla (key) con tolerancia fija */
    printf("\nQAG con epsrel = 1e-10 y distintas reglas Gauss-Kronrod\n");
    int keys[] = {GSL_INTEG_GAUSS15, GSL_INTEG_GAUSS31, GSL_INTEG_GAUSS61};
    const char *nk[] = {"GK15", "GK31", "GK61"};
    for (int k = 0; k < 3; k++) {
        double res, err;
        gsl_integration_qag(&F, a, b, 0.0, 1e-10, NMAX, keys[k], w, &res, &err);
        printf("  %-5s I = %.15f  err_est = %.2e  err_real = %.2e  subintervalos = %zu\n",
               nk[k], res, err, fabs(res - Iref), w->size);
    }

    /* QNG (no adaptativa) como contraste */
    {
        double res, err; size_t neval;
        int st = gsl_integration_qng(&F, a, b, 0.0, 1e-10, &res, &err, &neval);
        printf("  QNG   I = %.15f  err_est = %.2e  err_real = %.2e  nevals = %zu (%s)\n",
               res, err, fabs(res - Iref), neval, gsl_strerror(st));
    }
    gsl_integration_workspace_free(w);

    /*  EJERCICIO 2  */
    double x[200], y[200];
    int n = leer_csv("/home/carlos/MCEI_2620/Python/jupyter_lab/datos_sensor.csv", x, y, 200);
    if (n != 50) { fprintf(stderr, "Se esperaban 50 datos, hay %d\n", n); return 1; }
    double h = x[1] - x[0];
    for (int i = 1; i < n - 1; i++)
        if (fabs((x[i + 1] - x[i]) - h) > 1e-9) { fprintf(stderr, "No equiespaciado\n"); return 1; }
    printf("\n=== EJERCICIO 2 ===\n%d datos, h = %.4f, %d subintervalos\n", n, h, n - 1);

    double It = trapecio_datos(x, y, n);
    double Is = simpson13(y, 47, h) + simpson38(&y[46], h);   /* 46 + 3 intervalos */
    printf("Trapecio explícito           : %.10f\n", It);
    printf("Simpson 1/3 (46) + 3/8 (3)   : %.10f\n", Is);

    /* Representación interpolada: se construye un modelo continuo y se integra */
    const gsl_interp_type *tipos[] = {gsl_interp_linear, gsl_interp_cspline, gsl_interp_akima,
                                      gsl_interp_steffen};
    const char *nombres[] = {"lineal", "spline cúbico natural", "Akima", "Steffen"};
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    for (int k = 0; k < 4; k++) {
        gsl_spline *sp = gsl_spline_alloc(tipos[k], n);
        gsl_spline_init(sp, x, y, n);
        double Ii = gsl_spline_eval_integ(sp, x[0], x[n - 1], acc);
        printf("Interpolación %-22s: %.10f\n", nombres[k], Ii);
        gsl_spline_free(sp);
    }
    /*  EXTENSIÓN: derivadas  */
    printf("\n=== EXTENSIÓN: diferenciación ===\n");
    double dfw[200], dce[200];
    for (int i = 0; i < n - 1; i++) dfw[i] = (y[i + 1] - y[i]) / h;
    dfw[n - 1] = (y[n - 1] - y[n - 2]) / h;                   /* hacia atrás en el extremo */
    for (int i = 1; i < n - 1; i++) dce[i] = (y[i + 1] - y[i - 1]) / (2 * h);
    dce[0]     = (-3 * y[0] + 4 * y[1] - y[2]) / (2 * h);     /* unilateral O(h^2) */
    dce[n - 1] = (3 * y[n - 1] - 4 * y[n - 2] + y[n - 3]) / (2 * h);

    gsl_spline *cs = gsl_spline_alloc(gsl_interp_cspline, n);
    gsl_spline_init(cs, x, y, n);
    double emf = 0, emc = 0, ems = 0;
    FILE *out = fopen("derivadas_gsl.csv", "w");
    fprintf(out, "x,adelante,centrada,spline\n");
    for (int i = 0; i < n; i++) {
        double ex = 0.35 * 0.7 * cos(0.7 * x[i]) - 0.15 * 2.1 * sin(2.1 * x[i]) + 0.03;
        double ds = gsl_spline_eval_deriv(cs, x[i], acc);
        fprintf(out, "%.4f,%.8f,%.8f,%.8f\n", x[i], dfw[i], dce[i], ds);
        if (i < n - 1) emf = fmax(emf, fabs(dfw[i] - ex));
        emc = fmax(emc, fabs(dce[i] - ex));
        ems = fmax(ems, fabs(ds - ex));
    }
    fclose(out);
    printf("Error máx. adelante : %.3e\nError máx. centrada : %.3e\nError máx. spline   : %.3e\n",
           emf, emc, ems);
    gsl_spline_free(cs);
    gsl_interp_accel_free(acc);
    return 0;

}
