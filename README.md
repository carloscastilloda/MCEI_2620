Taller: integración y diferenciación numérica
Métodos computacionales en Ingeniería — MCEI_M  
Escuela Colombiana de Ingeniería Julio Garavito  
Elaborado por: Carlos Castillo

Comparación de tres entornos para integración y diferenciación numérica:
GNU Octave
C/C++ con GSL
Python con NumPy, SciPy y Matplotlib

El ejercicio 1 integra una función con primitiva conocida. El ejercicio 2 integra 50 datos equiespaciados. Una extensión analiza diferenciación numérica y sensibilidad al ruido.
---
Contenido
Estructura del repositorio
0Ejercicio 1: integración de una función
Ejercicio 2: integración de 50 datos equiespaciados
Extensión: diferenciación numérica
Comparación final
Reflexión
---
Estructura del repositorio
```text
taller/
├── octave/
│   └── taller_int_diff_CarlosCastillo.m          
├── C_C++/
    └integ_diff
│       ├── main_int_diff_CarlosCastillo.cpp             
│       └── CMakeLists.txt
└── python/
    └── jupiter_lab
        └── taller_integracion_diff_CarlosCastillo.ipynb         
        └── datos_sensor.csv              # 50 datos generados (Ejercicio 2)
```
Los tres programas leen `../datos_sensor.csv` a partir de la ruta /home/carlos/MCEI_2620/Python/jupyter_lab/datos_sensor.csv.


---
Ejercicio 1: integración de una función
A. Análisis matemático
1. Referencia analítica
La integral tiene primitiva cerrada. Separando:
$$
I=\int_0^8 e^{-0.4x},dx+0.5\int_0^8 e^{-0.4x}\sin(3x),dx
$$
Con
$$
\int_0^L e^{-ax},dx=\frac{1-e^{-aL}}{a},
\qquad
\int_0^L e^{-ax}\sin(bx),dx=\frac{b-e^{-aL}(a\sin bL+b\cos bL)}{a^2+b^2}
$$
y $a=0.4$, $b=3$, $L=8$:
$$
\boxed{I_{\mathrm{ref}} = 2.559824508302063}
$$
Verificado con `mpmath` a 30 dígitos: $2.55982450830206310809\ldots$


2. Interpretación del área acumulada
Si $f$ es una señal (potencia, caudal o concentración), $F(x)=\int_0^x f$ es la magnitud acumulada: energía entregada, volumen transportado o dosis total. Como $f>0$ en todo el dominio ($1+0.5\sin 3x\ge 0.5$), $F$ es monótona creciente y se satura porque la envolvente decae: el 76 % del área total se acumula antes de $x=3$.

3. Características que afectan la aproximación
Oscilación con periodo $2\pi/3\approx 2.09$: en $[0,8]$ hay $\sim 3.8$ ciclos. Con $n=10$ ($h=0.8$) hay menos de 3 puntos por ciclo, por lo que el trapecio submuestrea la oscilación.
La amplitud de la oscilación decae con la envolvente: el error se concentra en la región inicial, lo que favorece métodos adaptativos.
La función es $C^\infty$ y acotada. Se cumplen las hipótesis de suavidad: se espera orden $O(h^2)$ limpio para el trapecio y convergencia muy rápida de Gauss–Kronrod.

B. GNU Octave
Trapecio compuesto explícito (`trapecio_comp.m`). Tiempo promedio de 2000 repeticiones:
$n$	Integral aproximada	Error	Orden observado	Tiempo [s]
10	2.494072510185698	6.575e-02	—	4.4e-05
20	2.544929716211618	1.489e-02	2.14	4.4e-05
50	2.557502776878403	2.322e-03	2.03	3.9e-05
100	2.559246207984129	5.783e-04	2.01	4.0e-05
500	2.559801403459290	2.310e-05	2.00	4.8e-05
1000	2.559818732303387	5.776e-06	2.00	5.8e-05

El orden observado
$$
p=\log(e_k/e_{k+1})/\log(n_{k+1}/n_k)\to 2
$$
confirma $O(h^2)$. El tiempo casi no crece entre $n=10$ y $n=1000$ porque domina el overhead del intérprete (llamada a función anónima, `linspace`), no el número de evaluaciones.
Segunda aproximación: `integral` (adaptativa) y `quadgk`.
RelTol	`integral()`	Error real	Tiempo [s]
1e-6	2.559824508302063	4.4e-16	3.1e-04
1e-10	2.559824508302063	4.4e-16	3.1e-04
1e-14	2.559824508302063	0	6.0e-04
`quadgk` (RelTol $10^{-10}$): error estimado $6.4\times 10^{-12}$, error real $4.4\times 10^{-16}$.

C. C/C++ con GSL
Cuadratura adaptativa `gsl_integration_qag` (regla `GSL_INTEG_GAUSS21`, `epsabs = 0`, `limit = 1000`):
epsrel	Resultado	Error estimado	Error real	Subintervalos	Tiempo [s]
1e-4	2.559824508302063	1.458e-10	4.4e-16	2	1.3e-06
1e-6	2.559824508302063	1.458e-10	4.4e-16	2	1.3e-06
1e-10	2.559824508302063	1.458e-10	4.4e-16	2	1.3e-06
1e-13	2.559824508302063	2.842e-14	0	4	3.0e-06
Con `epsrel = 1e-10` y distintas reglas: GK15 usa 6 subintervalos; GK31 y GK61 resuelven con un solo intervalo. QNG (no adaptativa) converge con 43 evaluaciones.
Trapecio explícito en C: los mismos valores que Octave (idénticos hasta el último dígito), con tiempos de $2.6\times 10^{-7}$ s ($n=10$) a $1.6\times 10^{-5}$ s ($n=1000$). Aquí el tiempo sí escala linealmente con $n$, porque no hay overhead de intérprete.

Qué controla el programador y qué delega GSL
Controla el programador	Delega GSL
Definir `f(double x, void *params)` y empaquetar parámetros en un `struct`	Elección de nodos y pesos Gauss–Kronrod
Tolerancias `epsabs` / `epsrel` (criterio: $\lvert I-R\rvert\le\max(\mathrm{epsabs},\ \mathrm{epsrel}\lvert I\rvert)$)	

Estimación del error por diferencia Gauss vs. Kronrod
Regla (`key`: 15, 21, 31, 41, 51, 61 puntos)	Qué subintervalo bisecar (el de mayor error)
Reservar y liberar el `workspace`, fijar `limit`	Heurísticas de redondeo y detección de fallos
Revisar el código de estado (`GSL_SUCCESS`, `GSL_EMAXITER`, `GSL_EROUND`, …)	—
En el trapecio explícito el programador decide la malla y no obtiene estimación de error: hay que obtenerla por refinamiento (Richardson) o por comparación.

D. Python con SciPy
`scipy.integrate.quad` con `epsabs=0`:
epsrel	Resultado	Error estimado	Error real	neval	Tiempo [s]
1e-4	2.559824508302063	1.458e-10	4.4e-16	63	3.6e-05
1e-6	2.559824508302063	1.458e-10	4.4e-16	63	3.5e-05
1e-10	2.559824508302063	1.458e-10	4.4e-16	63	3.5e-05
1e-13	2.559824508302063	2.842e-14	0	147	8.1e-05

Valores por defecto de `quad`: `epsabs = epsrel = 1.49e-8`, `limit = 50`.


Comparación de interfaces SciPy vs. GSL
`quad` y `gsl_integration_qag` producen el mismo error estimado bit a bit ($1.458\times 10^{-10}$) y las mismas 63 evaluaciones ($= 21 \times 3$: el intervalo inicial y sus dos mitades), porque ambos descienden de QUADPACK (`quad` envuelve QAGS en Fortran; GSL es una reimplementación en C). La diferencia está en la interfaz, no en el algoritmo:
SciPy: una línea, cualquier callable, parámetros vía `args=`, sin gestión de memoria; devuelve `(valor, error)` y, con `full_output`, información diagnóstica. Los fallos se reportan como `IntegrationWarning`.
GSL: callback con firma fija, `workspace` explícito, elección de la regla y código de estado que el programador debe revisar.
Costo: GSL tarda $\sim 1.3,\mu\mathrm{s}$ frente a $\sim 35,\mu\mathrm{s}$ de SciPy, con el mismo número de evaluaciones. La diferencia ($\sim 27\times$) es el costo de llamar 63 veces a una función Python desde Fortran. Si $f$ fuera costosa (una simulación, por ejemplo), esa diferencia sería irrelevante.
---
Preguntas de análisis

1. ¿Qué cambia en el trapecio al aumentar $n$?  
El error disminuye como $h^2$: duplicar $n$ lo divide entre $\sim 4$, y multiplicar $n$ por 10 lo divide entre $\sim 100$ (de $5.8\times 10^{-4}$ con $n=100$ a $5.8\times 10^{-6}$ con $n=1000$). El costo crece linealmente ($n+1$ evaluaciones). Para ganar un dígito se necesita multiplicar $n$ por $\sqrt{10}\approx 3.2$; alcanzar $10^{-14}$ exigiría $n\sim 10^8$, y además el redondeo acumulado empezaría a dominar.
2. ¿Aumentar $n$ vs. reducir la tolerancia de un método adaptativo?  
Aumentar $n$ refina la malla uniformemente y a ciegas: gasta puntos por igual donde la función es suave y donde oscila, y no informa cuánto error queda. Reducir la tolerancia es fijar el objetivo (el error aceptable); el algoritmo decide dónde y cuánto refinar según su estimación local del error y se detiene al cumplirlo. La primera es una decisión sobre el procedimiento; la segunda, sobre el resultado.
3. ¿Mayor precisión implica necesariamente mayor costo?  
No. Entre métodos, el orden y la eficiencia pesan más que el número de puntos: QAG obtiene error de $4.4\times 10^{-16}$ con 63 evaluaciones, mientras el trapecio con 1001 evaluaciones se queda en $5.8\times 10^{-6}$ (diez órdenes de magnitud peor con 16 veces más evaluaciones). Dentro de un mismo método sí hay un compromiso: pasar de `epsrel` $10^{-10}$ a $10^{-13}$ subió el costo de 63 a 147 evaluaciones. Además, tolerancias por debajo de $\sim 10^{-15}$ (la precisión de `double`) no pueden cumplirse y solo generan advertencias de redondeo.
4. Error estimado vs. error real.  
Los estimadores son conservadores: $1.46\times 10^{-10}$ estimado frente a $4.4\times 10^{-16}$ real (un factor $\sim 10^5$ en GSL/SciPy) y $6.4\times 10^{-12}$ frente a $4.4\times 10^{-16}$ en `quadgk`. QUADPACK toma la diferencia entre las reglas Gauss (10 puntos) y Kronrod (21 puntos), que en realidad mide el error de la regla de Gauss —la menos precisa— y le aplica un escalado pesimista $\min(1,(200,\Delta)^{1.5})$. Esto explica por qué `epsrel` de $10^{-4}$, $10^{-6}$ y $10^{-10}$ dan el mismo resultado: el algoritmo ya “cree” cumplir las tres con 2 subintervalos. El error estimado es una cota práctica, no el error verdadero; es confiable para funciones suaves, pero puede subestimar el error si la función tiene singularidades o picos que los nodos no detectan.
---
Ejercicio 2: integración de 50 datos equiespaciados
A. Preparación
`datos_sensor.csv`: encabezado `x,y` y 50 filas. Los tres programas verifican que hay exactamente 50 datos, que no hay `NaN` y que $\max|\Delta x_i - h| < 10^{-9}$, con $h = 0.2$ e intervalo $[0,\ 9.8]$.
> **Observación clave:** 50 puntos equivalen a **49 subintervalos, un número impar**. Simpson 1/3 compuesto requiere un número par de subintervalos, así que **la condición de la regla no se cumple directamente**.

Opciones consideradas:
Simpson 1/3 en 46 intervalos + Simpson 3/8 en los últimos 3 (elegida): conserva $O(h^4)$ en todo el dominio.
Simpson 1/3 en 48 intervalos + trapecio en el último: más simple, pero el intervalo final degrada localmente el orden.
`scipy.integrate.simpson` (SciPy ≥ 1.11): aplica automáticamente la corrección de Cartwright en el último intervalo.
B–D. Resultados en los tres entornos
Método	Entorno	$I\approx\int_0^{9.8} y,dx$
Trapecio	Octave `trapz` / GSL explícito / SciPy `trapezoid`	21.1908463399
Simpson 1/3 (46) + 3/8 (3)	Octave / C / Python (manual)	21.1920398006
Simpson 1/3 (48) + trapecio (1)	Octave	21.1919529443
`simpson` (corrección Cartwright)	SciPy	21.1921130983
Interpolación lineal	GSL `gsl_interp_linear`	21.1908463399
Spline cúbico natural	GSL `cspline` / SciPy `CubicSpline`	21.1918882576
Akima	GSL	21.1920570526
Steffen (monótona)	GSL	21.1917706176
Cuando no se dispone de la función, el error del trapecio se puede estimar solo con los datos mediante Richardson: se compara $T_h$ con $T_{2h}$ (tomando uno de cada dos puntos). $(T_h-T_{2h})/3 \approx 1.1\times 10^{-3}$, que coincide en magnitud con la diferencia trapecio–Simpson ($1.19\times 10^{-3}$). Esto indica que el trapecio subestima la integral en $\sim 0.006,%$.
Validación a posteriori (solo para el informe; la función generadora no se usó en los cálculos): el valor exacto es $21.192018$. Errores: trapecio $1.2\times 10^{-3}$, Simpson 1/3+3/8 $2.2\times 10^{-5}$, Akima $3.9\times 10^{-5}$, SciPy `simpson` $9.5\times 10^{-5}$ y spline natural $1.3\times 10^{-4}$.

¿Son relevantes las diferencias?  
Numéricamente, Simpson es $\sim 50$ veces más preciso que el trapecio. En la práctica, la diferencia relativa ($0.006,%$) es mucho menor que la incertidumbre típica de un sensor real ($0.1$–$1,%$). Con estos datos, el trapecio es suficiente y el error de medición dominaría sobre el error de cuadratura. Simpson solo se justifica si los datos son muy precisos o si se necesita un orden alto por otras razones.
Integrar datos directamente vs. interpolarlos (GSL)
Integrar directamente (trapecio / Simpson)	Interpolar y luego integrar
Simple, rápido y transparente; no inventa información	Produce una función continua que se puede evaluar, derivar e integrar en cualquier subintervalo $[c,d]$
El trapecio es robusto ante ruido y datos no equiespaciados	Splines de orden alto logran precisión comparable a Simpson y admiten mallas no uniformes
Simpson exige malla uniforme y paridad de intervalos	Depende del modelo: condiciones de frontera (spline natural vs. not-a-knot), oscilaciones espurias y sobreajuste al ruido
Solo da la integral entre nodos	Más costo y más decisiones; puede sugerir una precisión que los datos no tienen
La interpolación lineal reproduce exactamente el trapecio: el trapecio es la integral de la interpolante lineal a trozos. De forma análoga, Simpson integra parábolas por tramos.
---
Extensión: diferenciación numérica
Se usaron diferencias hacia adelante, $O(h)$, y centradas, $O(h^2)$. Tratamiento de los extremos:
Hacia adelante: el último punto no tiene $x_{i+1}$; se deja como `NaN` (Octave, Python) o se usa diferencia hacia atrás (C).
Centrada: el primer y el último punto no tienen vecino a ambos lados. Se usan fórmulas unilaterales de segundo orden,
$$
f'_0\approx\frac{-3f_0+4f_1-f_2}{2h}
$$
y la simétrica en $x_n$ (`np.gradient(..., edge_order=2)` las implementa). El `gradient` por defecto de Octave usa en los extremos diferencias de primer orden.
Método	Error máximo (vs. derivada exacta, solo validación)
Adelante $O(h)$	8.07e-02
Centrada $O(h^2)$	1.67e-02
Derivada del spline cúbico (GSL `gsl_spline_eval_deriv`)	3.88e-02
El error del spline se concentra en los extremos, por la condición natural $y''=0$. Los tres entornos dan valores idénticos.
Sensibilidad al ruido
Se agregó ruido gaussiano $\sigma$ a los datos (promedio de 500 realizaciones, Python):
$\sigma$	Cambio relativo de la integral	Error relativo de la derivada (RMS / $\sigma$ de $y'$)
0.001	5.3e-05	1.6e-02
0.005	2.6e-04	7.7e-02
0.010	5.1e-04	1.6e-01
0.050	2.6e-03	7.8e-01
La derivada es $\sim 300$ veces más sensible que la integral.
Diferenciar amplifica el ruido. Para la diferencia centrada,
$$
\mathrm{Var}\left[\frac{\varepsilon_{i+1}-\varepsilon_{i-1}}{2h}\right]=\frac{\sigma^2}{2h^2}
$$
el error de ruido crece como $\sigma/h$. Al reducir $h$ baja el error de truncamiento $O(h^2)$, pero sube el de ruido. Existe un $h$ óptimo y refinar no siempre mejora el resultado.
Integrar promedia el ruido. Es una suma ponderada de muchos términos con errores independientes que se cancelan parcialmente. La varianza del trapecio es $\approx h^2 N\sigma^2$, así que la desviación estándar decrece como $\sigma\sqrt{hL}$ al refinar.
En frecuencia: derivar multiplica por $j\omega$ (filtro pasa-altos, realza el ruido de alta frecuencia); integrar divide por $j\omega$ (pasa-bajos, lo atenúa).
La figura `ext_python.png` muestra además que los extremos son los puntos más sensibles: las fórmulas unilaterales tienen coeficientes mayores ($-3$, $4$, $-1$), lo que amplifica más el ruido.
En la práctica, con señales reales conviene filtrar antes de derivar (Savitzky–Golay, suavizado por splines o un filtro pasa-bajos) o derivar un modelo ajustado en lugar de los datos crudos.
---
Comparación final
Criterio	Octave	C/C++ + GSL	Python + SciPy
Facilidad de implementación	Muy alta: vectorización nativa, `trapz` / `integral` en una línea	Baja: callbacks, `struct` de parámetros, memoria manual, lectura de CSV a mano	Muy alta: `quad`, `trapezoid`, `simpson` y `loadtxt`
Control del algoritmo	Medio: tolerancias, `Waypoints`, `MaxIntervalCount`	Máximo: regla GK, `limit`, `workspace`, códigos de estado, familia QNG/QAG/QAGS/CQUAD	Alto: tolerancias, `limit`, `points`, `weight`; la regla está fija
Manejo de tolerancias	`AbsTol` / `RelTol` por nombre	`epsabs` / `epsrel` explícitos; el programador revisa el `status`	`epsabs` / `epsrel` con valores por defecto de $1.49\times 10^{-8}$; avisos vía `warnings`
Integración de funciones	`integral`, `quadgk`, `quad`	QAG/QAGS/QAGI/CQUAD, Romberg, Gauss–Legendre fijo	`quad`, `fixed_quad`, `romb`, `quad_vec`
Integración de datos	`trapz`, `cumtrapz`; Simpson hay que programarlo	Trapecio manual; `gsl_spline_eval_integ` sobre interpolantes	`trapezoid`, `simpson` (maneja $n$ impar), `romb`, `CubicSpline.integrate`
Estimación del error	`quadgk` devuelve la cota	Devuelve la cota + código de estado	Devuelve la cota + `infodict` (neval, subintervalos)
Tiempo (Ej. 1)	Trapecio $n=1000$: $58,\mu\mathrm{s}$ · `integral`: $310,\mu\mathrm{s}$	Trapecio $n=1000$: $16,\mu\mathrm{s}$ · QAG: $1.3,\mu\mathrm{s}$	Trapecio $n=1000$: $34,\mu\mathrm{s}$ · `quad`: $35,\mu\mathrm{s}$
Visualización	Integrada, al estilo MATLAB	No tiene: se exporta CSV y se grafica fuera	Matplotlib, muy completa
Síntesis. Octave y Python sirven para prototipar y explorar. C + GSL da control y velocidad cuando la integral está dentro de un lazo costoso (optimización, Monte Carlo, simulación en tiempo real o un gemelo digital). En los tres casos el algoritmo de fondo (QUADPACK) es el mismo, así que la elección de entorno responde a necesidades de ingeniería de software, no de precisión.
---
Reflexión
1. De $f(x)$ a datos $(x_i,y_i)$
Con la función conocida se puede evaluar donde se quiera: el método elige los nodos (Gauss–Kronrod, adaptatividad), la tolerancia se puede fijar arbitrariamente y se obtienen estimaciones de error. Con datos se pierde todo eso: los nodos están fijos (la malla la decidió el experimento), no se puede refinar, no hay información entre muestras y aparece el ruido de medición. Por eso cambia la estrategia:
Se pasa de métodos adaptativos a reglas de Newton–Cotes compuestas (trapecio, Simpson) o a la integración de una interpolante.
Hay que verificar condiciones que antes eran automáticas: equiespaciado, paridad de intervalos (aquí 49, impar) y datos faltantes.
El error se estima con los mismos datos (Richardson con $2h$, comparación entre métodos), no con una referencia.
El método de orden alto deja de ser siempre preferible: si el ruido domina, un método robusto como el trapecio es tan bueno como Simpson, y el error de cuadratura es irrelevante frente a la incertidumbre del sensor.
2. Conceptos fundamentales
El flujo fue: (1) formular y obtener una referencia (analítica o de alta precisión); (2) prototipar en Octave para validar; (3) reimplementar en C/GSL y Python contra esa referencia; (4) comparar valor, error y tiempo.
Discretización. $h$ controla el error de truncamiento según el orden del método ($h^2$, $h^4$). Verificar el orden observado es la mejor prueba de que la implementación es correcta.
Error. Hay que distinguir error de truncamiento, de redondeo y de los datos. El error estimado por una biblioteca es una cota heurística, no el error real.
Tolerancia. Especifica el resultado deseado, no el procedimiento. Tiene un piso (la precisión de máquina) y conviene pedir solo lo que los datos o el modelo justifican.
Representación de datos. Una función, un vector de muestras y una interpolante son objetos distintos, con posibilidades distintas (evaluar, derivar, integrar en subintervalos). Las condiciones de cada regla (paridad, equiespaciado) dependen de esa representación.
Abstracción de bibliotecas. SciPy y GSL comparten algoritmo (QUADPACK) y difieren en cuánto exponen. Conocer lo que la abstracción oculta (la regla GK21, los valores por defecto de tolerancia, el `limit`) permite interpretar resultados sospechosos y ajustar el algoritmo cuando falla.
Validación cruzada. Obtener el mismo valor en tres entornos independientes da una confianza que ninguna implementación aislada ofrece.
