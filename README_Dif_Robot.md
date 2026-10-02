# Taller: Diferenciación numérica: velocidades de un robot diferencial
Métodos Computacionales en Ingeniería (MCEI_M)  
Escuela Colombiana de Ingeniería Julio Garavito  
Autor: Carlos Castillo
---

**Datos:** 51 muestras, h = 0.2 s, t ∈ [0, 10] s. Las trayectorias se generan desde las expresiones analíticas, lo que permite validar contra la solución exacta:

- ẋ = 0.16t + 0.18cos(0.45t), ẏ = 0.50 + 0.135 sin(0.45t)
- v = √(ẋ² + ẏ²), ω = (ẋÿ − ẏẍ)/(ẋ² + ẏ²) (curvatura × v)

---
Estructura del repositorio
```text
taller/
├── octave/
│   └── trayectoria_robot.csv 
│   └── taller_diff_robot_CarlosCastillo.m
│   └── graficas_robot.png
├── C_C++/
│   └── Diff_Robot/
│       ├── main_int_robot_CarlosCastillo.cpp
│       └── CMakeLists.txt
└── python/
    └── jupyter_lab/
        ├── taller_diff_robot_CarlosCastillo.ipynb
        └── graficas_python_robot.csv
        └── trayectoria_robot.csv
```
El archivo `trayectoria_robot.csv` fue generado por Octave

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

**Hallazgo importante:** el error de ω en i = 1 (5.9e-2) es 70 veces mayor que en el resto del interior. Se debe a que θ₀ se calcula con ẋ₀, ẏ₀ de primer orden (error O(h)); al hacer la diferencia centrada (θ₂ − θ₀)/2h, ese error O(h) se divide por 2h y queda un error O(1) en ω₁. Los extremos contaminan el primer punto interior en la segunda diferenciación encadenada. Lejos de los bordes el error de ω sí converge como O(h²) (8.4e-4 → 2.2e-4 → 5.4e-5 al dividir h por 2).

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

---

## Reflexión

**1. Reconstrucción del procedimiento.** (i) Se derivan x(t) e y(t) respecto a t por diferencias centrales → ẋᵢ, ẏᵢ (primera fuente de error: truncamiento O(h²) y amplificación de ruido ~σ/h; error mayor O(h) en extremos). (ii) Se combinan algebraicamente: vᵢ = √(ẋᵢ² + ẏᵢ²) y θᵢ = atan2(ẏᵢ, ẋᵢ) (el error se propaga; θ es mal condicionado si v → 0). (iii) Se desenvuelve θ (error si |ω|h ≥ π). (iv) Se deriva θ por diferencias centrales → ωᵢ (segunda diferenciación: el ruido se amplifica ~σ/h² y los errores de los extremos contaminan los puntos vecinos).

**2. Transferencia a datos reales.** Conservaría: la cadena (x,y) → (ẋ,ẏ) → (v,θ) → ω, el unwrap, la validación con una trayectoria sintética de solución conocida y el uso de diferencias sobre arreglos con paso posiblemente no uniforme (usar los tiempos reales, no un h nominal). Modificaría: filtrar o suavizar antes de derivar (Savitzky–Golay, spline de suavizado o filtro de Kalman con modelo cinemático), usar extremos de segundo orden, verificar timestamps duplicados o irregulares, y, si están disponibles, fusionar con odometría de ruedas (v = (v_R + v_L)/2, ω = (v_R − v_L)/b) o con el giróscopo de una IMU, que mide ω directamente. Criterio de selección: elegir el método y la ventana minimizando el error total estimado (truncamiento + ruido), a partir de la densidad espectral del ruido y la banda de interés de la dinámica (frecuencia de muestreo frente a frecuencia de corte), validando con datos sintéticos con el mismo nivel de ruido.

---

## Licencia

Material desarrollado con fines académicos para el curso **Métodos Computacionales en Ingeniería (MCEI_M)**, Escuela Colombiana de Ingeniería Julio Garavito.
