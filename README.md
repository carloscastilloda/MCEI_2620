# Taller: Integración y Diferenciación Numérica
Métodos Computacionales en Ingeniería (MCEI_M)  
Escuela Colombiana de Ingeniería Julio Garavito  
Autor: Carlos Castillo
---
## Descripción:

Este taller compara tres entornos computacionales para implementar métodos de integración y diferenciación numérica:
GNU Octave
C/C++ con GSL
Python con NumPy, SciPy y Matplotlib
Se abordan tres problemas:
Integración numérica de una función con primitiva conocida.
Integración de un conjunto de 50 datos equiespaciados.
Diferenciación numérica y análisis de sensibilidad al ruido.
---
Estructura del repositorio
```text
taller/
├── octave/
│   └── taller_int_diff_CarlosCastillo.m
├── C_C++/
│   └── integ_diff/
│       ├── main_int_diff_CarlosCastillo.cpp
│       └── CMakeLists.txt
└── python/
    └── jupyter_lab/
        ├── taller_integracion_diff_CarlosCastillo.ipynb
        └── datos_sensor.csv
```
El archivo `datos_sensor.csv` contiene los 50 datos generados para el Ejercicio 2.
> [!NOTE]
> Cada implementación usa la ruta relativa correspondiente para acceder al archivo de datos.
---
## Ejercicio 1: Integración de una función

### A. Análisis matemático

#### 1. Referencia analítica

La integral considerada posee primitiva cerrada:

$$
I = \int_0^8 e^{-0.4x}\,dx + 0.5\int_0^8 e^{-0.4x}\sin(3x)\,dx
$$

Se utilizan las integrales de referencia:

$$
\int_0^L e^{-ax}\,dx = \frac{1-e^{-aL}}{a}
$$

$$
\int_0^L e^{-ax}\sin(bx)\,dx = \frac{b-e^{-aL}\left[a\sin(bL)+b\cos(bL)\right]}{a^2+b^2}
$$

con $a=0.4$, $b=3$ y $L=8$, de donde se obtiene:

$$
\boxed{I_{\mathrm{ref}} = 2.559824508302063}
$$

El valor se verificó con `mpmath` usando 30 dígitos de precisión:

$$
I_{\mathrm{ref}} = 2.55982450830206310809\ldots
$$

#### 2. Interpretación del área acumulada

Si $f(x)$ representa una señal física (potencia, caudal o concentración), entonces

$$
F(x) = \int_0^x f(t)\,dt
$$

representa una magnitud acumulada: energía entregada, volumen transportado o dosis total.

Como $f(x) = e^{-0.4x}(1+0.5\sin 3x) > 0$, la función acumulada es monótona creciente. Además, por el decaimiento de la envolvente, aproximadamente el **76 % del área total se acumula antes de $x=3$**.

#### 3. Características que afectan la aproximación

- La oscilación tiene periodo $T = 2\pi/3 \approx 2.09$, por lo que en $[0,8]$ hay aproximadamente **3.8 ciclos**.
- Con $n=10$ ($h=0.8$) se obtienen menos de tres puntos por ciclo, de modo que el trapecio puede submuestrear la oscilación.
- La amplitud de la oscilación disminuye con la envolvente exponencial, por lo que el error tiende a concentrarse en la región inicial.
- La función es $C^\infty$ y acotada, por lo que se espera:
  - Trapecio compuesto: $O(h^2)$.
  - Gauss–Kronrod adaptativo: convergencia rápida.

### B. GNU Octave

Se utilizó el **trapecio compuesto explícito**. Resultados promediados sobre 2000 repeticiones:

|  $n$ | Integral aproximada |                Error | Orden observado |         Tiempo [s] |
| ---: | ------------------: | -------------------: | --------------: | -----------------: |
|   10 |   2.494072510185698 | $6.575\times10^{-2}$ |               — | $4.4\times10^{-5}$ |
|   20 |   2.544929716211618 | $1.489\times10^{-2}$ |            2.14 | $4.4\times10^{-5}$ |
|   50 |   2.557502776878403 | $2.322\times10^{-3}$ |            2.03 | $3.9\times10^{-5}$ |
|  100 |   2.559246207984129 | $5.783\times10^{-4}$ |            2.01 | $4.0\times10^{-5}$ |
|  500 |   2.559801403459290 | $2.310\times10^{-5}$ |            2.00 | $4.8\times10^{-5}$ |
| 1000 |   2.559818732303387 | $5.776\times10^{-6}$ |            2.00 | $5.8\times10^{-5}$ |

El orden observado se calcula como

$$
p = \frac{\log(e_k/e_{k+1})}{\log(n_{k+1}/n_k)}
$$

y tiende a $p \rightarrow 2$, lo que confirma el comportamiento $O(h^2)$.

#### Integración adaptativa

|     RelTol |      `integral()` |          Error real |         Tiempo [s] |
| ---------: | ----------------: | ------------------: | -----------------: |
|  $10^{-6}$ | 2.559824508302063 | $4.4\times10^{-16}$ | $3.1\times10^{-4}$ |
| $10^{-10}$ | 2.559824508302063 | $4.4\times10^{-16}$ | $3.1\times10^{-4}$ |
| $10^{-14}$ | 2.559824508302063 |         $\approx 0$ | $6.0\times10^{-4}$ |

`quadgk` con `RelTol = 1e-10` reportó un error estimado de $6.4\times10^{-12}$, frente a un error real de $4.4\times10^{-16}$: la estimación es conservadora.

### C. C/C++ con GSL

Se utilizó `gsl_integration_qag` con la regla `GSL_INTEG_GAUSS21`, `epsabs = 0` y `limit = 1000`.

|   `epsrel` |         Resultado |        Error estimado |          Error real | Subintervalos |         Tiempo [s] |
| ---------: | ----------------: | --------------------: | ------------------: | ------------: | -----------------: |
|  $10^{-4}$ | 2.559824508302063 | $1.458\times10^{-10}$ | $4.4\times10^{-16}$ |             2 | $1.3\times10^{-6}$ |
|  $10^{-6}$ | 2.559824508302063 | $1.458\times10^{-10}$ | $4.4\times10^{-16}$ |             2 | $1.3\times10^{-6}$ |
| $10^{-10}$ | 2.559824508302063 | $1.458\times10^{-10}$ | $4.4\times10^{-16}$ |             2 | $1.3\times10^{-6}$ |
| $10^{-13}$ | 2.559824508302063 | $2.842\times10^{-14}$ |         $\approx 0$ |             4 | $3.0\times10^{-6}$ |

#### ¿Qué controla el programador y qué delega GSL?

| Control del programador             | Responsabilidad de GSL                  |
| ----------------------------------- | --------------------------------------- |
| Definir `f(double x, void *params)` | Elección de nodos y pesos Gauss–Kronrod |
| Definir `epsabs` y `epsrel`         | Estimación del error                    |
| Seleccionar la regla                | Selección del subintervalo a refinar    |
| Reservar/liberar el `workspace`     | Heurísticas de redondeo                 |
| Fijar `limit`                       | Detección de fallos                     |
| Revisar los códigos de estado       | —                                       |

En el trapecio explícito, en cambio, el programador controla completamente la malla pero no obtiene directamente una estimación del error.

### D. Python con SciPy

Se utilizó `scipy.integrate.quad` con `epsabs=0`.

|   `epsrel` |         Resultado |        Error estimado |          Error real | Evaluaciones |         Tiempo [s] |
| ---------: | ----------------: | --------------------: | ------------------: | -----------: | -----------------: |
|  $10^{-4}$ | 2.559824508302063 | $1.458\times10^{-10}$ | $4.4\times10^{-16}$ |           63 | $3.6\times10^{-5}$ |
|  $10^{-6}$ | 2.559824508302063 | $1.458\times10^{-10}$ | $4.4\times10^{-16}$ |           63 | $3.5\times10^{-5}$ |
| $10^{-10}$ | 2.559824508302063 | $1.458\times10^{-10}$ | $4.4\times10^{-16}$ |           63 | $3.5\times10^{-5}$ |
| $10^{-13}$ | 2.559824508302063 | $2.842\times10^{-14}$ |         $\approx 0$ |          147 | $8.1\times10^{-5}$ |

Valores por defecto de `quad`:

$$ epsabs = epsrel = 1.49e-8 $$

limit  = 50


### E. Comparación SciPy vs. GSL

`quad` y `gsl_integration_qag` producen el mismo error estimado ($1.458\times10^{-10}$) y usan las mismas **63 evaluaciones** en este caso. La diferencia principal está en la interfaz:

| SciPy                                      | GSL                                                      |
| ------------------------------------------ | -------------------------------------------------------- |
| Interfaz sencilla                          | Interfaz de bajo nivel                                   |
| Acepta cualquier `callable`                | Callback con firma fija                                  |
| Parámetros mediante `args=`                | Parámetros mediante `struct`                             |
| Sin gestión manual de memoria              | `workspace` explícito                                    |
| Devuelve `(valor, error)`                  | Devuelve resultado + código de estado                    |
| Advertencias con `IntegrationWarning`      | Códigos `GSL_SUCCESS`, `GSL_EMAXITER`, `GSL_EROUND`      |

En tiempo de ejecución, $t_{\mathrm{GSL}} \approx 1.3\ \mu s$ frente a $t_{\mathrm{SciPy}} \approx 35\ \mu s$: una diferencia de aproximadamente **27×**, atribuible principalmente al costo de las llamadas entre Python y Fortran.

---

## Ejercicio 2: Integración de 50 datos equiespaciados

### A. Preparación

El archivo `datos_sensor.csv` contiene:

- Encabezado: `x,y`
- 50 observaciones
- Paso $h = 0.2$
- Intervalo $[0,\ 9.8]$

Verificaciones realizadas:

```text
n = 50
NaN = 0
max|Δx_i - h| < 1e-9
```

### B–D. Resultados

| Método                 | Entorno              | $I \approx \int_0^{9.8} y\,dx$ |
| ---------------------- | -------------------- | -----------------------------: |
| Trapecio               | Octave / GSL / SciPy |              **21.1908463399** |
| Simpson 1/3 + 3/8      | Octave / C / Python  |              **21.1920398006** |
| Simpson 1/3 + trapecio | Octave               |              **21.1919529443** |
| `scipy.integrate.simpson` | SciPy             |              **21.1921130983** |
| Interpolación lineal   | GSL                  |              **21.1908463399** |
| Spline cúbico natural  | GSL / SciPy          |              **21.1918882576** |
| Akima                  | GSL                  |              **21.1920570526** |
| Steffen                | GSL                  |              **21.1917706176** |

Valor exacto usado para la validación: $I_{\mathrm{exacto}} = 21.192018$

#### Errores

| Método                    |              Error |
| ------------------------- | -----------------: |
| Trapecio                  | $1.2\times10^{-3}$ |
| Simpson 1/3 + 3/8         | $2.2\times10^{-5}$ |
| Akima                     | $3.9\times10^{-5}$ |
| `scipy.integrate.simpson` | $9.5\times10^{-5}$ |
| Spline natural            | $1.3\times10^{-4}$ |

### Integración directa vs. interpolación

| Integración directa                         | Interpolación + integración                 |
| ------------------------------------------- | ------------------------------------------- |
| Simple y rápida                             | Produce una función continua                |
| No introduce información adicional          | Permite evaluar en cualquier subintervalo   |
| Robusta ante ruido                          | Permite derivar e integrar                  |
| El trapecio funciona con mallas no uniformes | Puede trabajar con mallas no uniformes     |
| Simpson requiere condiciones específicas    | Depende del modelo de interpolación         |
| Menor costo computacional                   | Mayor costo y número de decisiones          |

> La interpolación lineal reproduce exactamente la regla del trapecio, porque el trapecio equivale a integrar la interpolante lineal a trozos.

---

## Extensión: Diferenciación numérica

Esquemas utilizados:

- Diferencias hacia adelante: $O(h)$.
- Diferencias centradas: $O(h^2)$.
- En los extremos, fórmulas unilaterales de segundo orden. Por ejemplo, en el primer punto:

$$
f'_0 \approx \frac{-3f_0 + 4f_1 - f_2}{2h}
$$

| Método                     |        Error máximo |
| -------------------------- | ------------------: |
| Adelante $O(h)$            | $8.07\times10^{-2}$ |
| Centrada $O(h^2)$          | $1.67\times10^{-2}$ |
| Derivada del spline cúbico | $3.88\times10^{-2}$ |

El error del spline se concentra en los extremos, debido a la condición de frontera natural $y''=0$.

### Sensibilidad al ruido

Se agregó ruido gaussiano con desviación estándar $\sigma$ y se realizaron **500 realizaciones**.

| $\sigma$ | Cambio relativo de la integral | Error relativo de la derivada |
| -------: | -----------------------------: | ----------------------------: |
|    0.001 |             $5.3\times10^{-5}$ |            $1.6\times10^{-2}$ |
|    0.005 |             $2.6\times10^{-4}$ |            $7.7\times10^{-2}$ |
|    0.010 |             $5.1\times10^{-4}$ |            $1.6\times10^{-1}$ |
|    0.050 |             $2.6\times10^{-3}$ |            $7.8\times10^{-1}$ |

En este experimento, la derivada resulta aproximadamente **300 veces más sensible al ruido** que la integral.

Para la diferencia centrada:

$$
\mathrm{Var}\left[
\frac{\varepsilon_{i+1}-\varepsilon_{i-1}}{2h}
\right] =
\frac{\sigma^{2}}{2h^{2}}
$$

por lo que el error asociado al ruido crece aproximadamente como $\sigma/h$. Esto genera un compromiso entre el error de truncamiento, el error por ruido y el tamaño de paso $h$: **reducir $h$ no siempre mejora la estimación de la derivada**.

### Interpretación en frecuencia

- La diferenciación multiplica el espectro por $j\omega$: actúa como un **filtro pasa-altos** y amplifica las componentes de alta frecuencia.
- La integración divide por $j\omega$: se comporta como un **filtro pasa-bajos**.

---

## Reflexión

### De $f(x)$ a datos $(x_i, y_i)$

Cuando se conoce la función $f(x)$, es posible evaluarla en cualquier punto, y los métodos adaptativos pueden elegir los nodos y ajustar la tolerancia.

Con datos experimentales, en cambio:

- Los nodos los determina el experimento.
- No es posible refinar la malla.
- No hay información entre muestras.
- Aparece el ruido de medición.

La estrategia cambia entonces hacia el trapecio, Simpson o la interpolación, y el error se estima comparando métodos o usando el refinamiento disponible en los datos.

---

## Licencia

Material desarrollado con fines académicos para el curso **Métodos Computacionales en Ingeniería (MCEI_M)**, Escuela Colombiana de Ingeniería Julio Garavito.
