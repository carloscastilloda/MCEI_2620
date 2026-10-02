/* =====================================================================
 * Taller Diferenciación Robot
 * MCEI_M - Escuela Colombiana de Ingeniería
 * Elaborador Por: Carlos Castillo
 * ===================================================================== */

 /*diferencias centrales sobre datos tabulados
 * y demostración conceptual de gsl_deriv_central sobre una función evaluable.*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics_double.h>

#define NMAX 10000

/* Derivada por diferencias: centrada en el interior, unilateral en los extremos */
static void diff_central(const double *t, const double *f, double *df, int n) {
    for (int i = 1; i < n - 1; i++)
        df[i] = (f[i + 1] - f[i - 1]) / (t[i + 1] - t[i - 1]);
    df[0]     = (f[1] - f[0]) / (t[1] - t[0]);
    df[n - 1] = (f[n - 1] - f[n - 2]) / (t[n - 1] - t[n - 2]);
}

/* Desenvolvimiento angular: elimina saltos mayores a pi */
static void unwrap(double *th, int n) {
    double offset = 0.0;
    for (int i = 1; i < n; i++) {
        double d = th[i] + offset - th[i - 1];
        if (d > M_PI)       offset -= 2.0 * M_PI * ceil((d - M_PI) / (2.0 * M_PI));
        else if (d < -M_PI) offset += 2.0 * M_PI * ceil((-d - M_PI) / (2.0 * M_PI));
        th[i] += offset;
    }
}

/* Funciones evaluables para gsl_deriv_central */
static double x_fun(double t, void *p) { (void)p; return 0.08 * t * t + 0.40 * sin(0.45 * t); }
static double y_fun(double t, void *p) { (void)p; return 0.50 * t + 0.30 * (1.0 - cos(0.45 * t)); }

int main(int argc, char **argv) {
    const char *fname = argc > 1 ? argv[1] : "trayectoria_robot.csv";
    FILE *fp = fopen(fname, "r");
    if (!fp) { perror(fname); return 1; }

    static double t[NMAX], x[NMAX], y[NMAX], vx[NMAX], vy[NMAX], v[NMAX], th[NMAX], w[NMAX];
    char line[256];
    int n = 0;
    if (!fgets(line, sizeof line, fp)) { fclose(fp); return 1; }        /* encabezado */
    while (n < NMAX && fscanf(fp, "%lf,%lf,%lf", &t[n], &x[n], &y[n]) == 3) n++;
    fclose(fp);

    struct timespec a, b;
    clock_gettime(CLOCK_MONOTONIC, &a);
    diff_central(t, x, vx, n);
    diff_central(t, y, vy, n);
    for (int i = 0; i < n; i++) {
        v[i]  = gsl_hypot(vx[i], vy[i]);
        th[i] = atan2(vy[i], vx[i]);
    }
    unwrap(th, n);
    diff_central(t, th, w, n);
    clock_gettime(CLOCK_MONOTONIC, &b);
    double t_exec = (b.tv_sec - a.tv_sec) + 1e-9 * (b.tv_nsec - a.tv_nsec);

    /* Validación contra la solución analítica */
    double ev_int = 0, ew_int = 0, ev = 0, ew = 0;
    for (int i = 0; i < n; i++) {
        double T = t[i];
        double xd = 0.16 * T + 0.18 * cos(0.45 * T), yd = 0.50 + 0.135 * sin(0.45 * T);
        double xdd = 0.16 - 0.081 * sin(0.45 * T),   ydd = 0.06075 * cos(0.45 * T);
        double dv = fabs(v[i] - hypot(xd, yd));
        double dw = fabs(w[i] - (xd * ydd - yd * xdd) / (xd * xd + yd * yd));
        if (dv > ev) ev = dv;
        if (dw > ew) ew = dw;
        if (i > 0 && i < n - 1) { if (dv > ev_int) ev_int = dv; if (dw > ew_int) ew_int = dw; }
    }
    printf("C/GSL  | tiempo: %.3e s  (n = %d)\n", t_exec, n);
    printf("Error max interior  v: %.3e   omega: %.3e\n", ev_int, ew_int);
    printf("Error max global    v: %.3e   omega: %.3e\n", ev, ew);
    printf("Media v = %.6f m/s, desv. est. = %.6f (gsl_stats)\n",
           gsl_stats_mean(v, 1, n), gsl_stats_sd(v, 1, n));

    /* gsl_deriv_central: requiere una FUNCIÓN evaluable en puntos arbitrarios t +/- h/2, t +/- h */
    gsl_function Fx = {&x_fun, NULL}, Fy = {&y_fun, NULL};
    double T = 5.0, dx, dy, ex, ey;
    gsl_deriv_central(&Fx, T, 1e-3, &dx, &ex);
    gsl_deriv_central(&Fy, T, 1e-3, &dy, &ey);
    printf("\ngsl_deriv_central en t = %.1f: xdot = %.10f (+/- %.1e), ydot = %.10f (+/- %.1e)\n",
           T, dx, ex, dy, ey);
    printf("Exacto:                       xdot = %.10f,            ydot = %.10f\n",
           0.16 * T + 0.18 * cos(0.45 * T), 0.50 + 0.135 * sin(0.45 * T));
    printf("Tabulado (dif. centrada h=0.2): xdot = %.10f,            ydot = %.10f\n", vx[25], vy[25]);

    FILE *fo = fopen("resultados_c.csv", "w");
    for (int i = 0; i < n; i++)
        fprintf(fo, "%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n", t[i], vx[i], vy[i], v[i], th[i], w[i]);
    fclose(fo);
    return 0;
}
