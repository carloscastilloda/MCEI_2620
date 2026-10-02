<<<<<<< HEAD
# Taller: Integración y Diferenciación Numérica
=======
Taller: Diferenciación numérica: velocidades de un robot diferencial
>>>>>>> 1adc783 (Taller diferenciación numérica: robot diferencial en Octave, Python y C++/GSL)
Métodos Computacionales en Ingeniería (MCEI_M)  
Escuela Colombiana de Ingeniería Julio Garavito  
Autor: Carlos Castillo
---
<<<<<<< HEAD
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
=======
**Datos:** 51 muestras, h = 0.2 s, t ∈ [0, 10] s. Las trayectorias se generan desde las expresiones analíticas, lo que permite validar contra la solución exacta:

- ẋ = 0.16t + 0.18cos(0.45t), ẏ = 0.50 + 0.135 sin(0.45t)
- v = √(ẋ² + ẏ²), ω = (ẋÿ − ẏẍ)/(ẋ² + ẏ²) (curvatura × v)

---
## Archivos

| Archivo | Contenido |
|---|---|
| `trayectoria_robot.csv` | Archivo generado desde Octave (t,x,y) |
| `taller_diff_robot_CarlosCastillo.m` | Diferencias centrales explícitas, `unwrap`, gráficas: `graficas_robot.png` |
| `taller_diff_robot_CarlosCastillo.ipyn` | Versión vectorizada con `np.gradient`, comparación con Octave y estudio de h con ruido → `graficas_python_robot.png` |
| `main_diff_robot_CarlosCastillo.cpp`, `CMakeLists.txt` | Diferencias centrales en C, unwrap propio, `gsl_stats`, demostración de `gsl_deriv_central` |
---
## Resultados

| Métrica | Octave | Python | C/GSL |
|---|---|---|---|
| Error máx. v, interior | 2.34e-4 | 2.34e-4 | 2.34e-4 |
| Error máx. ω, interior (i = 1…N−2) | 5.90e-2 | 5.90e-2 | 5.90e-2 |
| Error máx. ω, i = 2…N−3 | 8.4e-4 | 8.4e-4 | 8.4e-4 |
| Error máx. v / ω, incluyendo extremos | 2.29e-2 / 1.36e-1 | ídem | ídem |
| Tiempo de cálculo (51 muestras) | ~1.1 ms | ~0.6 ms | ~24 µs |

Los tres entornos coinciden hasta 1e-10 (diferencia atribuible solo al redondeo de escritura del CSV). Con `np.gradient(..., edge_order=2)` el error global de ω baja de 1.36e-1 a 1.07e-2.

El error de ω en i = 1 (5.9e-2) es 70 veces mayor que en el resto del interior. Se debe a que θ₀ se calcula con ẋ₀, ẏ₀ de primer orden (error O(h)); al hacer la diferencia centrada (θ₂ − θ₀)/2h, ese error O(h) se divide por 2h y queda un error O(1) en ω₁. Los extremos contaminan el primer punto interior en la segunda diferenciación encadenada. Lejos de los bordes el error de ω sí converge como O(h²) (8.4e-4 a 2.2e-4 a 5.4e-5 al dividir h por 2).

---
## Tabla comparativa

| Criterio | Octave | Python (NumPy) | C/C++ + GSL |
|---|---|---|---|
| Carga de datos | `dlmread` en una línea | `np.loadtxt` en una línea | `fopen`/`fscanf` manual, gestión de memoria y encabezado |
| Cálculo de ẋ, ẏ | Diferencias explícitas con indexación vectorial | `np.gradient` (centrada interior, unilateral en extremos) | Bucle explícito sobre arreglos |
| Cálculo de v | `sqrt(vx.^2+vy.^2)` | `np.sqrt`/`np.hypot` | `gsl_hypot` en bucle |
| Cálculo de θ | `unwrap(atan2(...))` incorporado | `np.unwrap(np.arctan2(...))` | `atan2` + unwrap implementado a mano |
| Cálculo de ω | Diferencias explícitas sobre θ | `np.gradient(theta, t)` | Misma función de diferencias reutilizada |
| Manejo de arreglos | Vectores nativos, índices desde 1 | `ndarray`, broadcasting, índices desde 0 | Arreglos estáticos C, tamaño fijo |
| Facilidad de implementación | Alta | Muy alta | Media-baja |
| Control sobre el algoritmo | Alto (fórmula visible) | Medio (extremos ocultos en `gradient`, controlables con `edge_order`) | Total |
| Tiempo de ejecución | ~1 ms (intérprete) | ~0.6 ms | ~24 µs (compilado, -O2) |

**Qué reemplaza NumPy:** `np.gradient` reemplaza los bucles/indexaciones de diferencias explícitas. En el interior usa diferencias centrales (con fórmula para paso no uniforme si `t` lo es); en los extremos usa por defecto diferencias unilaterales de primer orden (`edge_order=1`), idénticas a las de Octave en este taller; con `edge_order=2` usa fórmulas unilaterales de tres puntos O(h²).

**Sobre `gsl_deriv_central`:** recibe un `gsl_function` y evalúa f en t±h/2 y t±h, aplicando extrapolación de Richardson y estimando el error. En t = 5 s da ẋ = 0.6869287479 (error estimado 8e-11), exacto a 10 cifras, mientras la diferencia centrada tabulada con h = 0.2 da 0.68708. Con datos tabulados no existe una función que evaluar en puntos arbitrarios: solo hay muestras en la rejilla fija, así que la estrategia correcta es operar sobre los arreglos (o interpolar/ajustar un spline y derivar esa función, a costa de introducir el modelo de interpolación).

---

## Preguntas de análisis

**1. v desde ẋ, ẏ vs. desde diferencias de posición.** Calcular v = √(ẋ² + ẏ²) con derivadas centradas estima la rapidez *instantánea* en tᵢ con error O(h²). Calcular directamente ‖pᵢ₊₁ − pᵢ‖/h da la longitud de la cuerda dividida por h: es la rapidez *media* en el intervalo, asociada a tᵢ + h/2 y con error O(h) si se atribuye a tᵢ; además la cuerda siempre subestima la longitud de arco en tramos curvos. La primera opción también entrega la dirección (necesaria para θ y ω); la segunda solo la magnitud.

**2. Por qué un pequeño error en ẋ, ẏ afecta más a ω.** ω es una segunda derivada disfrazada: ω = (ẋÿ − ẏẍ)/v². Numéricamente se diferencian dos veces los datos, y cada diferenciación divide el error por un factor ~h. Además, el error angular es δθ ≈ δv⊥/v, de modo que cuando v es pequeño (arranques, giros cerrados) un error fijo en ẋ, ẏ produce un error grande en θ. El taller lo muestra: error relativo en v ~0.02 %, en ω hasta ~1 % en el interior y mucho mayor junto a los extremos.

**3. Reducir h con ruido fijo.** El error total tiene dos términos: truncamiento ~C·h² y ruido amplificado ~σ/h (para v) y ~σ/h² (para ω). Con σ = 1 mm (RMS, puntos interiores):

| h (s) | err v sin ruido | err v con ruido | err ω sin ruido | err ω con ruido |
|---|---|---|---|---|
| 0.4 | 6.8e-4 | 1.7e-3 | 1.2e-3 | 5.1e-3 |
| 0.2 | 1.7e-4 | 4.3e-3 | 3.2e-4 | 2.0e-2 |
| 0.1 | 4.2e-5 | 6.3e-3 | 8.4e-5 | 8.0e-2 |
| 0.05 | 1.0e-5 | 1.5e-2 | 2.1e-5 | 2.6e-1 |
| 0.01 | 4.2e-7 | 7.0e-2 | 8.7e-7 | 6.7 |

Sin ruido, reducir h mejora como h²; con ruido, por debajo de un h óptimo el error crece y ω se vuelve inutilizable (6.7 rad/s de error con h = 0.01 s). Existe un h óptimo, aproximadamente h* ~ (σ/C)^(1/3) para v.

**4. Desenvolvimiento de θ.** `atan2` devuelve valores en (−π, π]. Si el robot cruza la dirección θ = ±π, aparece un salto artificial de 2π entre muestras consecutivas y la diferencia (θᵢ₊₁ − θᵢ₋₁)/2h produce un pico espurio de ~2π/2h (≈ 15.7 rad/s con h = 0.2). `unwrap` suma múltiplos de 2π cuando el salto supera π, restaurando la continuidad. En esta trayectoria θ ∈ [0.23, 1.20] rad y no hay saltos, pero el paso es imprescindible en general (p. ej. trayectorias circulares). Requiere que el robot gire menos de π por muestra, es decir |ω|·h < π.

---

## Conclusiones

1. **Sobre el método:** la diferencia centrada da O(h²) en el interior y el resultado es idéntico en los tres entornos; la exactitud depende del esquema numérico, del tratamiento de los extremos y del encadenamiento de derivadas, no de la herramienta. Los extremos de primer orden degradan no solo los bordes sino el primer punto interior de ω, por lo que conviene usar fórmulas unilaterales de segundo orden o descartar esos puntos.
2. **Sobre el entorno:** Octave y Python permiten expresar el algoritmo en pocas líneas vectorizadas y son ideales para prototipar y graficar; C/GSL es ~40 veces más rápido y da control total (útil en un controlador embebido en tiempo real), a cambio de implementar a mano carga de datos, unwrap y manejo de memoria. Las rutinas `gsl_deriv_*` están pensadas para funciones evaluables, no para muestras.

## Reflexión

**1. Reconstrucción del procedimiento.** (i) Se derivan x(t) e y(t) respecto a t por diferencias centrales → ẋᵢ, ẏᵢ (primera fuente de error: truncamiento O(h²) y amplificación de ruido ~σ/h; error mayor O(h) en extremos). (ii) Se combinan algebraicamente: vᵢ = √(ẋᵢ² + ẏᵢ²) y θᵢ = atan2(ẏᵢ, ẋᵢ) (el error se propaga; θ es mal condicionado si v → 0). (iii) Se desenvuelve θ (error si |ω|h ≥ π). (iv) Se deriva θ por diferencias centrales → ωᵢ (segunda diferenciación: el ruido se amplifica ~σ/h² y los errores de los extremos contaminan los puntos vecinos).

**2. Transferencia a datos reales.** Conservaría: la cadena (x,y) → (ẋ,ẏ) → (v,θ) → ω, el unwrap, la validación con una trayectoria sintética de solución conocida y el uso de diferencias sobre arreglos con paso posiblemente no uniforme (usar los tiempos reales, no un h nominal). Modificaría: filtrar o suavizar antes de derivar (Savitzky–Golay, spline de suavizado o filtro de Kalman con modelo cinemático), usar extremos de segundo orden, verificar timestamps duplicados o irregulares, y, si están disponibles, fusionar con odometría de ruedas (v = (v_R + v_L)/2, ω = (v_R − v_L)/b) o con el giróscopo de una IMU, que mide ω directamente. Criterio de selección: elegir el método y la ventana minimizando el error total estimado (truncamiento + ruido), a partir de la densidad espectral del ruido y la banda de interés de la dinámica (frecuencia de muestreo frente a frecuencia de corte), validando con datos sintéticos con el mismo nivel de ruido.

>>>>>>> 1adc783 (Taller diferenciación numérica: robot diferencial en Octave, Python y C++/GSL)
