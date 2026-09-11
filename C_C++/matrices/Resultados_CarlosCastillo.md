# Taller de Álgebra Lineal Numérica — Métodos Computacionales en Ingeniería (MCEI)

**Universidad Escuela Colombiana de Ingeniería Julio Garavito**
**Docente:** Ing. Alexander Pérez Ruiz MSc. PhD.
**Estudiante:** Carlos Alberto Castillo Daza

Implementación dual en **GNU Octave** y **C++ con Eigen3**.

# PARTE 1 — Problemas propuestos

## Punto 1. Inversa de Moore-Penrose

$$A=\begin{bmatrix}1&0&2\\2&-1&5\\0&1&-1\\1&3&-1\end{bmatrix}\in\mathbb{R}^{4\times3}$$

### 1.1 A es deficiente de rango

La tercera columna es combinación lineal de las otras dos:

$$c_3 = 2c_1 - c_2 \qquad\Rightarrow\qquad \operatorname{rank}(A)=2 < 3$$

Los valores singulares lo confirman:

| $\sigma_i$ | valor |
|---|---|
| $\sigma_1$ | 6.134996470821 |
| $\sigma_2$ | 3.218977835123 |
| $\sigma_3$ | **5.067452e-16** ← cero numérico |

Tolerancia de rango $\max(m,n)\cdot\sigma_1\cdot\varepsilon = 5.449\times10^{-15}$, luego
$\operatorname{rank}_{\text{num}}(A)=2$.

**Consecuencia:** $A^TA$ es singular ($\det(A^TA) = -2.31\times10^{-13}$,
$\operatorname{cond}_2(A^TA)=6.27\times10^{16}$) y la fórmula
$A^+=(A^TA)^{-1}A^T$ **no aplica**. Si se fuerza igual, el resultado es basura:

$$\|(A^TA)^{-1}A^T - A^+\|_F = 3.24\times10^{-1}$$

es decir, un error del orden de la propia magnitud de la respuesta. La única vía
correcta es la SVD con truncación:

$$A = U\Sigma V^T \quad\Rightarrow\quad A^+ = V\Sigma^+U^T,\qquad
\Sigma^+_{ii}=\begin{cases}1/\sigma_i & \sigma_i>\text{tol}\\ 0 & \text{en otro caso}\end{cases}$$

### 1.2 Resultado

$$A^{+}=\begin{bmatrix}
0.05128205 & 0.07692308 & 0.02564103 & 0.12820513\\
0.04358974 & 0.01538462 & 0.07179487 & 0.25897436\\
0.05897436 & 0.13846154 & -0.02051282 & -0.00256410
\end{bmatrix}$$

Coincide con `pinv(A)` de Octave y con
`A.completeOrthogonalDecomposition().pseudoInverse()` de Eigen
($\|\cdot\|_F = 1.49\times10^{-16}$).

### 1.3 Verificación de las cuatro propiedades de Penrose

| # | Propiedad | Residuo $\|\cdot\|_F$ | Veredicto |
|---|---|---|---|
| 1 | $AA^{+}A = A$ | 2.433e-15 | CUMPLE |
| 2 | $A^{+}AA^{+} = A^{+}$ | 1.428e-16 | CUMPLE |
| 3 | $(AA^{+})^{T}=AA^{+}$ | 2.642e-16 | CUMPLE |
| 4 | $(A^{+}A)^{T}=A^{+}A$ | 4.228e-16 | CUMPLE |

Las cuatro se cumplen al nivel del épsilon de máquina ($\varepsilon=2.22\times10^{-16}$).

### 1.4 Interpretación geométrica

$AA^{+}$ y $A^{+}A$ son **proyectores ortogonales** (simétricos e idempotentes)
sobre $\operatorname{Im}(A)$ y $\operatorname{Im}(A^T)$ respectivamente. Su traza
es igual al rango:

$$\operatorname{tr}(AA^{+}) = \operatorname{tr}(A^{+}A) = 2 = \operatorname{rank}(A)$$

$$A^{+}A=\begin{bmatrix}0.3333&0.3333&0.3333\\0.3333&0.8333&-0.1667\\0.3333&-0.1667&0.8333\end{bmatrix}\ne I_3$$

Es decir: por ser deficiente de rango, $A^+$ no es inversa ni por izquierda ni
por derecha. Sigue siendo la matriz que da la solución de mínimos cuadrados de
norma mínima, que es lo que garantizan las cuatro propiedades.

---
## Punto 2. Sistemas de ecuaciones lineales (10×10)

El enunciado entrega $A$ pero no $b$. Se construye
$b = A\cdot\mathbf{1}$ para conocer la solución exacta ($x=\mathbf{1}$) y poder
medir el **error** y no solamente el residual:

$$b^T = [16,\;16,\;18,\;18,\;24,\;19,\;16,\;20,\;22,\;20]$$

### 2.1 Caracterización previa de A

| Propiedad | Valor |
|---|---|
| $\det(A)$ | 358 816 |
| $\operatorname{cond}_2(A)$ | 51.7150 |
| $\operatorname{cond}_\infty(A)$ | 89.7448 |
| ¿Simétrica? | Sí |
| Autovalores extremos | $-3.1471 \ldots 19.3216$ |

$A$ es simétrica pero indefinida (hay autovalores negativos), así que
Cholesky no es aplicable . Con $\operatorname{cond}_2=51.7$
la matriz está muy bien condicionada.

### 2.2 Resultados de los tres métodos

Los tres recuperan $x=\mathbf{1}$ con 8 cifras exactas.

| Método | $\|Ax-b\|_2$ | error relativo | tiempo medio (20 000 rep.) | costo |
|---|---|---|---|---|
| Gauss con pivoteo parcial (a mano) | 5.024e-15 | 6.101e-16 | 7.07e-07 s | $\tfrac{2}{3}n^3$ |
| LU (`PartialPivLU`) | 5.024e-15 | 4.399e-16 | 8.04e-07 s | $\tfrac{2}{3}n^3$ |
| QR (`HouseholderQR`) | 1.137e-14 | 1.013e-15 | 2.59e-06 s | $\tfrac{4}{3}n^3$ |
| `FullPivLU` (referencia) | 5.617e-15 | 1.530e-15 | — | $\tfrac{2}{3}n^3$ + pivoteo |

Cota teórica: $\text{err}_{\text{rel}}\le \operatorname{cond}_2(A)\cdot\varepsilon = 1.148\times10^{-14}$.
Los cuatro métodos quedan un orden de magnitud por debajo de la cota.

### 2.3 Calidad de las factorizaciones

| Comprobación | Valor |
|---|---|
| $\|PA-LU\|_F$ | 2.026e-15 |
| $\|A-QR\|_F$ | 8.156e-15 |
| $\|Q^TQ-I\|_F$ | 1.147e-15 |
| Permutaciones de fila en Gauss | 5 |

`diag(U)` $= [4, 3.5, 4.1429, 4.5517, 3.0682, 7.8543, -1.9345, 2.6719, -5.0992, -2.1399]$ —
ningún pivote pequeño, lo que confirma el buen condicionamiento.

### 2.4 Comparación de estabilidad y eficiencia

- **Estabilidad.** Los tres son *backward stable*. QR es el más robusto en
  general porque $Q$ es ortogonal y no amplifica el error
  ($\operatorname{cond}_2(Q)=1$, verificado: $\|Q^TQ-I\|_F\approx10^{-15}$),
  mientras que LU depende del factor de crecimiento de los pivotes. Aquí, con
  $\operatorname{cond}=52$, la diferencia es irrelevante.
- **Eficiencia.** QR cuesta el doble de operaciones que LU
  ($\frac{4}{3}n^3$ vs $\frac{2}{3}n^3$) y eso se refleja en el tiempo medido:
  2.59 µs vs 0.80 µs, un factor ≈3.2.
- Múltiples lados derechos. La ventaja decisiva de LU: se factoriza una
  vez en $O(n^3)$ y cada nuevo $b$ cuesta solo $O(n^2)$. El código lo
  demuestra resolviendo un segundo $b$ reutilizando la misma factorización.

**Conclusión:** para este sistema, LU con pivoteo parcial es la elección
correcta (más barata, igual de precisa). QR solo se justifica si la matriz
estuviera mal condicionada o el sistema fuera rectangular.

---
## Punto 3. Condicionamiento numérico y estabilidad

Matriz mal condicionada de referencia: **Hilbert**, $H_{ij}=1/(i+j-1)$.

### 3.1 Número de condición

| $n$ | $\operatorname{cond}_2(H_n)$ | dígitos perdidos | dígitos útiles (de 16) |
|---|---|---|---|
| 2 | 1.928e+01 | 1.3 | 14.7 |
| 4 | 1.551e+04 | 4.2 | 11.8 |
| 6 | 1.495e+07 | 7.2 | 8.8 |
| 8 | 1.526e+10 | 10.2 | 5.8 |
| **10** | **1.603e+13** | **13.2** | **2.8** |
| 12 | 1.747e+16 | 16.2 | **-0.2** |

Regla práctica: se pierden $\log_{10}(\operatorname{cond})$ cifras significativas.
Para $n=12$ no queda ninguna cifra correcta en doble precisión.

Valores singulares de $H_{10}$: desde $1.752$ hasta $1.093\times10^{-13}$ — los
últimos están ya al nivel del ruido de redondeo.

### 3.2 Sistema original vs. perturbado

Perturbación determinista $\delta b = 10^{-10}\mathbf{1}$:

| Cantidad | Valor |
|---|---|
| Perturbación relativa en $b$ | 6.714e-11 |
| Cambio relativo en $x$ | 5.676e-04 |
| Factor de amplificación real | **8.45e+06** |
| Cota teórica $\operatorname{cond}_2(H)$ | 1.603e+13 |

Una perturbación de 11 dígitos por debajo de $b$ produce un cambio en la
cuarta cifra decimal de $x$: 7 órdenes de magnitud de amplificación. La cota
teórica se cumple holgadamente (es una cota de peor caso).

Además, incluso sin perturbar hay error: el solo hecho de almacenar $H$ en
doble precisión ya es una perturbación de $\sim10^{-16}$, y el error relativo de
$x_0$ frente a $\mathbf{1}$ es **4.12e-05**.

**Contraste decisivo** — la *misma* perturbación aplicada a la matriz bien
condicionada del punto 2:

| Matriz | $\operatorname{cond}_2$ | pert. rel. en $b$ | cambio rel. en $x$ | amplificación |
|---|---|---|---|---|
| Hilbert(10) | 1.60e+13 | 6.71e-11 | 5.68e-04 | **8.45e+06** |
| $A$ del punto 2 | 5.17e+01 | 5.24e-12 | 2.29e-11 | **4.36** |

### 3.3 Estabilidad de los métodos sobre el problema mal condicionado

| Método | error rel. en $x$ | $\|Hx-b\|$ | sensibilidad a la perturbación |
|---|---|---|---|
| LU (Gauss) | 2.160e-04 | 2.937e-16 | 3.682e-04 |
| QR | 3.385e-04 | 2.937e-16 | 2.569e-04 |
| SVD completa | 4.496e-03 | 2.488e-15 | 6.428e-04 |
| **SVD truncada (rango 9/10)** | **3.134e-06** | 2.147e-15 | — |

### 3.4 Discusión

1. **Condicionamiento ≠ estabilidad.** El condicionamiento es una propiedad del
   problema; la estabilidad, del algoritmo. Ningún algoritmo estable
   puede arreglar un problema mal condicionado.
2. **Residual pequeño ≠ solución correcta.** Los residuales son $\sim10^{-16}$ en
   los tres métodos (son *backward stable*: resuelven exactamente un problema
   vecino) y sin embargo el error en $x$ es $10^{-4}$. Usar el residual como
   único criterio de calidad es un error clásico.
3. **La SVD truncada gana por dos órdenes de magnitud**: al descartar
   $\sigma_{10}=1.09\times10^{-13}$ se elimina la dirección que amplifica el
   ruido. Es regularización (Tikhonov / TSVD) y es la respuesta práctica al mal
   condicionamiento.

---
## Punto 4. Inversa por método directo, QR y SVD

Se prueba sobre dos matrices para que la comparación sea informativa.

### Caso 1 — $A$ del punto 2 (bien condicionada, $\operatorname{cond}_2=51.7$)

| Método | $\|AA^{-1}-I\|_F$ | $\|A^{-1}A-I\|_F$ | dif. rel. vs `inv()` | tiempo |
|---|---|---|---|---|
| Directo (`inverse`) | 2.828e-15 | 6.658e-15 | — | 7.15e-05 s |
| QR | 4.218e-15 | 4.591e-15 | 8.00e-16 | 1.39e-05 s |
| SVD | 8.291e-15 | 8.432e-15 | 1.66e-15 | 3.80e-05 s |

Error al resolver $Ax=b$ con cada inversa: **2.2e-15 – 4.2e-15** en los tres.
Cota teórica $\operatorname{cond}\cdot\varepsilon=1.15\times10^{-14}$.
Los tres métodos son equivalentes; decide el costo.

### Caso 2 — Hilbert(10) (mal condicionada, $\operatorname{cond}_2=1.60\times10^{13}$)

| Método | $\|AA^{-1}-I\|_F$ | $\|A^{-1}A-I\|_F$ | error al resolver $Ax=b$ |
|---|---|---|---|
| Directo (`inverse`) | 1.684e-04 | 6.754e-03 | 6.589e-03 |
| QR | 3.397e-04 | 3.644e-03 | 3.606e-03 |
| SVD | 7.546e-03 | 7.809e-03 | 4.839e-03 |
| **SVD truncada (rango 8/10)** | 1.414e+00 **(!)** | 1.414e+00 **(!)** | **5.455e-06** |

Cota teórica: $\operatorname{cond}\cdot\varepsilon = 3.56\times10^{-3}$.

Para la pseudoinversa truncada, $AA^{+}$ es un proyector de rango $k$, así que 
$\|AA^{+}-I\|_F=\sqrt{n-k}=\sqrt2=1.414$ por construcción, 
no por imprecisión. La métrica que importa es el error en la
solución: ahí la SVD truncada es 1000 veces mejor que cualquier otra.

### Conclusiones

1. **Bien condicionada:** los tres coinciden a $\sim10^{-15}$; gana la inversa
   directa por costo ($2n^3$ frente a $\sim\frac{4}{3}n^3+n^3$ de QR y
   $12\text{–}22n^3$ de la SVD).
2. **Mal condicionada:**  los tres errores quedan
   del orden de $\operatorname{cond}(A)\cdot\varepsilon$, que es exactamente la
   cota teórica. Esa pérdida es del problema, no del algoritmo.
3. La ventaja real de la SVD no es un error más pequeño en la inversa
   completa, sino que es el único método que diagnostica (muestra los
   $\sigma$ pequeños) y permite regularizar truncándolos.
4. Para resolver $Ax=b$ nunca se calcula $A^{-1}$.
   `x = A\b` es más barato ($\frac{2}{3}n^3$ vs $2n^3+n^2$) y más preciso.
   La inversa explícita solo se justifica cuando se necesita la matriz en sí,
   como la matriz de covarianza $(X^TX)^{-1}\sigma^2$ en estadística.
---
## Punto 5. Consulta a un LLM sobre mínimos cuadrados

### 5.1 Prompt construido con la plantilla del enunciado

> **Actúa como** profesor universitario de ingeniería biomédica y estudiante de
> doctorado en ingeniería con énfasis en modelos computacionales, con
> experiencia en procesamiento de señales biomédicas e implementación de
> algoritmos matriciales en Octave y C++/Eigen.
>
> **El contexto es** un laboratorio de análisis de movimiento. Se registra
> electromiografía de superficie (sEMG) del bíceps y el tríceps braquial junto
> con el torque isométrico de flexión de codo medido en un dinamómetro. A partir
> de ventanas de 250 ms se extraen características de la señal —valor RMS de
> cada músculo, valor absoluto medio (MAV) y frecuencia mediana del espectro— y
> se desea estimar el torque con un modelo lineal multivariable
> $\tau = \beta_0+\beta_1\,\text{RMS}_{bi}+\beta_2\,\text{RMS}_{tri}+\beta_3\,\text{MAV}_{bi}+\beta_4\,f_{med}$
> ajustado sobre $n$ ventanas. Las características extraídas de un mismo canal
> son fuertemente colineales, de modo que la matriz de diseño está mal
> condicionada y los estudiantes deben entender no solo cómo obtener
> $\hat\beta$ sino por qué unos algoritmos son numéricamente preferibles a
> otros y qué hacer cuando el problema mismo está mal planteado.
>
> **Necesito que desarrolles** una explicación completa del método de mínimos
> cuadrados para el ajuste de modelos lineales: planteamiento del problema de
> optimización, deducción de las ecuaciones normales, interpretación geométrica
> como proyección ortogonal, y las tres vías de solución (ecuaciones normales,
> factorización QR y SVD/pseudoinversa) con su costo computacional y su
> comportamiento frente al número de condición.
>
> **El resultado debe estar en formato** documento técnico estructurado con
> secciones numeradas, ecuaciones en LaTeX y una tabla comparativa de los tres
> métodos.
>
> **Incluye** los supuestos del teorema de Gauss-Markov y cuáles de ellos viola
> típicamente una señal biomédica, la relación
> $\operatorname{cond}(X^TX)=\operatorname{cond}(X)^2$ y sus implicaciones
> prácticas —ilustrada con el ajuste polinómico para remoción de deriva de línea
> de base en ECG—, las métricas de bondad de ajuste (residual, MSE, RMSE, $R^2$)
> junto con la validación cruzada, y los criterios para decidir qué método usar
> en cada situación.

### 5.2 Respuesta

#### 5.2.1 Planteamiento

Con $n$ ventanas y $p$ parámetros ($n>p$), el sistema $X\beta=y$ es
**sobredeterminado**: en general no tiene solución exacta porque
$y\notin\operatorname{Im}(X)$ —el torque medido nunca es una combinación lineal
exacta de las características de sEMG, y no debería serlo: el sEMG es una señal
estocástica y el dinamómetro tiene su propio ruido. Se busca entonces el
$\hat\beta$ que minimice la norma euclídea del residual:

$$\hat\beta=\arg\min_{\beta\in\mathbb{R}^p}\;\|y-X\beta\|_2^2$$

Se usa la norma 2 (y no la 1 o la $\infty$) por tres razones: es diferenciable,
produce un problema lineal en $\beta$, y bajo ruido gaussiano coincide con el
estimador de máxima verosimilitud. En señales biomédicas hay una cuarta razón
práctica: la norma 2 admite formulación recursiva (RLS), que es lo que permite
correr el estimador en línea en una prótesis mioeléctrica.

#### 5.2.2 Deducción de las ecuaciones normales

Sea $S(\beta)=\|y-X\beta\|_2^2=(y-X\beta)^T(y-X\beta)=y^Ty-2\beta^TX^Ty+\beta^TX^TX\beta$.

$$\frac{\partial S}{\partial\beta}=-2X^Ty+2X^TX\beta=0
\quad\Longrightarrow\quad \boxed{X^TX\hat\beta=X^Ty}$$

La hessiana $\partial^2S/\partial\beta^2=2X^TX$ es semidefinida positiva, y
definida positiva si $X$ tiene rango completo por columnas; el punto crítico es
por tanto un **mínimo global único**.

#### 5.2.3 Interpretación geométrica

La condición de optimalidad equivale a $X^T(y-X\hat\beta)=X^Tr=0$: el residual
es **ortogonal al espacio columna de $X$**. Es decir, $\hat y=X\hat\beta$ es la
**proyección ortogonal** de $y$ sobre $\operatorname{Im}(X)$:

$$\hat y = X(X^TX)^{-1}X^Ty = Hy,\qquad H=H^T=H^2 \;(\text{matriz sombrero})$$

Esta es la comprobación numérica más útil del ajuste: si $\|X^Tr\|_\infty$ no es
prácticamente cero, el $\hat\beta$ calculado está mal. La lectura de ingeniería
es directa: el modelo se queda con la parte del torque que las características
de sEMG **pueden** explicar, y descarta el resto como residual.

#### 5.2.4 Las tres vías de solución

**(a) Ecuaciones normales.** Se forma $X^TX$ ($p\times p$) y se resuelve con
Cholesky (posible porque $X^TX$ es simétrica definida positiva).
Costo $\approx np^2+\frac{1}{3}p^3$ — el más barato.
**Problema fatal:** $\operatorname{cond}_2(X^TX)=\operatorname{cond}_2(X)^2$.
Se **duplica** el número de dígitos perdidos. Con $\operatorname{cond}(X)\ge10^8$
el resultado no tiene ninguna cifra correcta.

**(b) Factorización QR.** $X=QR$ con $Q^TQ=I_p$ y $R$ triangular superior.
Como la norma 2 es invariante bajo transformaciones ortogonales,

$$\|y-X\beta\|_2^2=\|Q^Ty-R\beta\|_2^2+\|(I-QQ^T)y\|_2^2$$

y el mínimo se alcanza resolviendo el sistema triangular $R\hat\beta=Q^Ty$.
Costo $\approx 2np^2-\frac{2}{3}p^3$ (≈2× las normales).
**Ventaja:** nunca se forma $X^TX$, así que el condicionamiento efectivo es
$\operatorname{cond}(X)$ y no su cuadrado. Es el método **por defecto** (lo que
hace `A\b` en Octave/MATLAB).

**(c) SVD / pseudoinversa.** $X=U\Sigma V^T$, $\hat\beta=V\Sigma^+U^Ty=X^+y$.
Costo $\approx 4np^2+8p^3$ — el más caro.
**Ventaja:** funciona aunque $X$ sea deficiente de rango (devuelve la solución
de norma mínima), expone el espectro completo y permite regularizar por
truncación o Tikhonov. Es la vía indicada cuando las características de sEMG son
redundantes, que es casi siempre.

| Criterio | Ec. normales | QR | SVD |
|---|---|---|---|
| Costo (flops) | $np^2+\tfrac13p^3$ | $2np^2-\tfrac23p^3$ | $4np^2+8p^3$ |
| Condicionamiento efectivo | $\operatorname{cond}(X)^2$ | $\operatorname{cond}(X)$ | $\operatorname{cond}(X)$ |
| ¿Soporta rango deficiente? | No | Con pivoteo de columnas | **Sí** |
| ¿Permite regularizar? | Ridge | Limitado | **Sí (TSVD)** |
| ¿Diagnostica el problema? | No | Parcialmente | **Sí** |
| Cuándo usarlo | $p$ pequeño y $X$ bien condicionada | **Por defecto** | Rango dudoso, colinealidad |

#### 5.2.5 Supuestos de Gauss-Markov y cuáles viola una señal biomédica

Si (i) $\mathbb{E}[\varepsilon]=0$, (ii) $\operatorname{Var}(\varepsilon)=\sigma^2I$
(homocedasticidad y no correlación), y (iii) $X$ es de rango completo y no
aleatoria, entonces $\hat\beta$ es el **BLUE** (*Best Linear Unbiased
Estimator*): insesgado y de mínima varianza entre los estimadores lineales
insesgados, con

$$\operatorname{Var}(\hat\beta)=\sigma^2(X^TX)^{-1},\qquad
\hat\sigma^2=\frac{\text{SSE}}{n-p}$$

**No** se requiere normalidad para Gauss-Markov; la normalidad se necesita solo
para las pruebas $t$ y $F$ y los intervalos de confianza.

En sEMG los tres supuestos se violan de formas concretas y conocidas, y conviene
nombrarlas antes de reportar un $p$-valor:

1. **Heterocedasticidad.** La varianza del sEMG crece con el nivel de
   activación: a mayor torque, mayor amplitud *y* mayor dispersión. El supuesto
   $\operatorname{Var}(\varepsilon)=\sigma^2I$ falla, y la corrección es mínimos
   cuadrados ponderados (WLS) o una transformación estabilizadora de varianza.
2. **Autocorrelación.** Si las ventanas de 250 ms se calculan con solapamiento
   (lo habitual, 50–75 %), las filas de $X$ **no son independientes** y los
   errores estándar salen artificialmente pequeños. La solución es GLS con
   blanqueo ($\Omega^{-1/2}$), o evaluar con ventanas disjuntas.
3. **Errores en las variables.** $X$ no está libre de ruido: el RMS es él mismo
   un estimador con varianza. Esto sesga $\hat\beta$ hacia cero (*atenuación de
   regresión*), y es la justificación formal de usar ventanas suficientemente
   largas: reducen la varianza de la característica a costa de resolución
   temporal.
4. **Rango completo.** Es el supuesto que falla de manera más evidente, y se
   trata en la sección siguiente.

#### 5.2.6 $\operatorname{cond}(X^TX)=\operatorname{cond}(X)^2$ en la práctica

**Caso 1: características redundantes de sEMG.** Para una señal de media cero y
distribución aproximadamente gaussiana, MAV y RMS están ligados analíticamente:

$$\mathbb{E}|x| = \sigma\sqrt{2/\pi},\qquad \text{RMS}=\sigma
\quad\Longrightarrow\quad \frac{\text{MAV}}{\text{RMS}}=\sqrt{\tfrac{2}{\pi}}\approx0.7979$$

Es decir, MAV y RMS del mismo canal son **casi proporcionales**: su correlación
empírica en sEMG real suele estar por encima de 0.98. Incluir ambas en $X$ crea
colinealidad severa, $\sigma_{\min}(X)\to0$ y $\operatorname{cond}(X)$ explota.
Los $\hat\beta$ individuales se vuelven inestables (grandes, de signos opuestos,
cambiando de sujeto a sujeto) aunque la **predicción** $\hat y$ siga siendo
buena — el síntoma clásico de la colinealidad.

**Caso 2: deriva de línea de base en ECG por ajuste polinómico.** Ajustar un
polinomio de grado $g$ en el tiempo para sustraer la deriva es un problema de
mínimos cuadrados con matriz de **Vandermonde**, que es el ejemplo canónico de
mal condicionamiento. Para un registro de 10 s a 1 kHz ($n=10\,000$):

| Grado $g$ | $\operatorname{cond}(X)$, base cruda $t^k$ | $\operatorname{cond}(X^TX)$ | $\operatorname{cond}(X)$ con $t$ normalizado a $[-1,1]$ |
|---|---|---|---|
| 3 | 1.96e+03 | 3.86e+06 | **8.22** |
| 5 | 4.60e+05 | 2.12e+11 | **43.2** |
| 8 | 2.10e+09 | **4.41e+18** | **553** |

La lectura es contundente. Con grado 8 y base cruda, las ecuaciones normales
trabajan con $\operatorname{cond}(X^TX)=4.4\times10^{18} > 1/\varepsilon$: el
sistema es **numéricamente singular** y no queda ninguna cifra correcta, aunque
el problema matemático esté perfectamente bien planteado. La misma QR sobre la
misma base ya conserva 7 dígitos. Y el simple cambio de variable
$t\mapsto 2(t-t_{\min})/(t_{\max}-t_{\min})-1$ baja el condicionamiento de
$2\times10^9$ a **553**: siete órdenes de magnitud **gratis**, sin cambiar el
modelo ni un solo dato. Usar una base ortogonal (Legendre, Chebyshev) lo lleva
a $\operatorname{cond}=1$.

**Ejemplo mínimo de por qué formar $X^TX$ destruye información** (Läuchli): con
$\epsilon=10^{-9}$,
$X=\begin{bmatrix}1&1\\ \epsilon&0\\0&\epsilon\end{bmatrix}$ tiene rango 2, pero
$X^TX=\begin{bmatrix}1+\epsilon^2&1\\1&1+\epsilon^2\end{bmatrix}$ se redondea a
$\begin{bmatrix}1&1\\1&1\end{bmatrix}$ porque $1+\epsilon^2=1$ en doble
precisión: la matriz se vuelve singular y la información del rango se pierde de
forma irreversible.

#### 5.2.7 Métricas de bondad de ajuste y validación

$$r=y-X\hat\beta,\quad \text{SSE}=\|r\|_2^2,\quad \text{MSE}=\frac{\text{SSE}}{n},\quad
\text{RMSE}=\sqrt{\text{MSE}}$$
$$R^2=1-\frac{\text{SSE}}{\text{SST}},\qquad
R^2_{\text{adj}}=1-(1-R^2)\frac{n-1}{n-p}$$

El RMSE está en las **unidades físicas** de $y$ (N·m en este caso) y es la
métrica interpretable: un RMSE de 2 N·m sobre torques de 0–40 N·m significa un
error del 5 % de la escala. $R^2$ es adimensional y siempre crece al añadir
regresores, por lo que para comparar modelos de distinto tamaño debe usarse
$R^2_{\text{adj}}$. En estadística se prefiere el **estimador insesgado**
$\hat\sigma^2=\text{SSE}/(n-p)$ sobre $\text{SSE}/n$.

**Advertencia específica de modelos biomédicos:** estas métricas calculadas
sobre los mismos datos del ajuste miden **memorización**, no capacidad
predictiva. En un modelo que se va a usar sobre un sujeto nuevo, la métrica
válida es el error en validación cruzada, y la partición debe ser **por sujeto**
(*leave-one-subject-out*), nunca por ventana: ventanas contiguas del mismo
sujeto son casi idénticas, y mezclarlas entre entrenamiento y prueba infla el
$R^2$ de manera espuria. Es el error metodológico más frecuente en la literatura
de sEMG.

#### 5.2.8 Criterio de decisión práctico

1. Calcular $\operatorname{cond}_2(X)$ **antes** de elegir el método. Es una
   línea de código y evita horas de depuración.
2. $\operatorname{cond}(X)<10^4$ → cualquier método sirve; QR por seguridad.
3. $10^4<\operatorname{cond}(X)<10^8$ → QR obligatorio; las ecuaciones normales
   ya pierden la mitad de las cifras.
4. $\operatorname{cond}(X)>10^8$ → SVD, y **antes de eso**, revisar el
   planteamiento: centrar y escalar las columnas, normalizar el eje temporal,
   usar una base ortogonal, o eliminar características redundantes (no incluir
   RMS y MAV del mismo canal). Casi siempre el mal condicionamiento es un
   síntoma de un modelo mal parametrizado, no un problema de aritmética.
5. Si la colinealidad es intrínseca al fenómeno y no se puede eliminar,
   regularizar: TSVD, *ridge* ($\hat\beta_\lambda=(X^TX+\lambda I)^{-1}X^Ty$) o
   reducción previa por PCA sobre las características.
6. Siempre: verificar $\|X^Tr\|_\infty\approx0$ al final y reportar el error en
   validación cruzada, no el del ajuste.

---
# PARTE 2 — Problema práctico: convertidor DC–DC

Datos: `datos_convertidor_realista.csv`, $n=100$ mediciones en régimen estacionario.

| Variable | media | desv. est. | mín | máx |
|---|---|---|---|---|
| Voltaje $V$ [V] | 21.4032 | 5.9497 | 12.11 | 31.74 |
| Corriente $I$ [A] | 6.9785 | 2.9311 | 2.07 | 11.86 |
| Temperatura $T$ [°C] | 50.8799 | 14.6712 | 25.25 | 74.50 |
| Potencia $P$ [W] | 134.2667 | 68.5366 | 24.36 | 332.96 |

## Entregable 1 — Coeficientes del modelo

$X=[\mathbf{1}\;V\;I\;T]\in\mathbb{R}^{100\times4}$, $y=P$.

| Coeficiente | Ec. normales | QR | SVD | ColPivQR |
|---|---|---|---|---|
| $\beta_0$ (intercepto) | -123.41834970 | -123.41834970 | -123.41834970 | -123.41834970 |
| $\beta_1$ ($V$) | 6.32256854 | 6.32256854 | 6.32256854 | 6.32256854 |
| $\beta_2$ ($I$) | 19.04206847 | 19.04206847 | 19.04206847 | 19.04206847 |
| $\beta_3$ ($T$) | -0.20682478 | -0.20682478 | -0.20682478 | -0.20682478 |

$$\boxed{\;\hat P = -123.4183 + 6.3226\,V + 19.0421\,I - 0.2068\,T\;}$$

Diferencias respecto a QR: ecuaciones normales **3.70e-12**, con `inv()`
explícita **2.47e-12**, SVD **1.43e-14**, ColPivQR **9.97e-14**.

Los cuatro métodos coinciden en las 8 cifras mostradas porque
$\operatorname{cond}(X)=347$ es moderado. **Pero la jerarquía ya se ve**: las
ecuaciones normales tienen un error 2–3 órdenes de magnitud mayor que QR/SVD.
Con $\operatorname{cond}(X)\sim10^7$ esas mismas ecuaciones normales no darían
ninguna cifra correcta.

## Entregable 2 — Influencia relativa de cada variable

| Variable | $\beta$ | error est. | $t$ | $\beta$ estandarizado | influencia |
|---|---|---|---|---|---|
| $V$ | 6.322569 | 0.267704 | **23.62** | 0.548865 | **39.0 %** |
| $I$ | 19.042068 | 0.548927 | **34.69** | 0.814361 | **57.9 %** |
| $T$ | -0.206825 | 0.109681 | **-1.89** | -0.044274 | **3.1 %** |
| intercepto | -123.418350 | 9.590755 | -12.87 | — | — |

Los $\beta$ crudos **no son comparables** porque cada variable tiene unidades y
escalas distintas. La comparación correcta usa coeficientes estandarizados
$\beta_j^{*}=\beta_j\,s_{x_j}/s_y$, que responden a "¿cuántas desviaciones
estándar cambia $P$ por cada desviación estándar de $x_j$?".

Efecto de recorrer **todo el rango** medido de cada variable:

| Variable | rango | efecto sobre $P$ |
|---|---|---|
| $V$ | 19.63 V | **+124.11 W** |
| $I$ | 9.79 A | **+186.42 W** |
| $T$ | 49.25 °C | **-10.19 W** |

Correlaciones con $P$: $V$ **0.5228**, $I$ **0.8021**, $T$ **-0.1841**.
Correlaciones entre regresores: $V$–$I$ -0.034, $V$–$T$ -0.038, $I$–$T$ -0.146
→ **prácticamente ortogonales, no hay colinealidad** (diseño experimental
factorial bien hecho).

## Entregable 3 — Residual y error cuadrático medio

| Métrica | Valor |
|---|---|
| $\|r\|_2$ | 155.0408 W |
| SSE | 24 037.6517 W² |
| **MSE** | **240.3765 W²** |
| **RMSE** | **15.5041 W** |
| MAE | 12.4479 W |
| $\max\|r\|$ | 47.9571 W |
| $R^2$ | **0.94830942** |
| $R^2$ ajustado | 0.94669409 |
| RMSE / media($P$) | 11.55 % |
| media($r$) | -2.98e-14 **[OK]** |
| $\|X^Tr\|_\infty$ | 1.52e-10 **[OK]** |

Las dos últimas filas son las **verificaciones numéricas** del ajuste: el
residual tiene media nula y es ortogonal a las columnas de $X$, tal como exige
la condición de optimalidad de mínimos cuadrados.

El modelo explica el **94.8 %** de la varianza, pero un RMSE de 15.5 W sobre una
media de 134 W (11.5 %) es alto para un convertidor de potencia: hay
**estructura no capturada**, no solo ruido. La sección de diagnóstico lo aclara.

## Entregable 4 — Comparación Octave ↔ C++/Eigen3

| Cantidad | C++ / Eigen3 | Referencia LAPACK (Octave) | Diferencia |
|---|---|---|---|
| $\beta_0$ | -123.41834970 | -123.41834970 | < 1e-9 |
| $\beta_1$ | 6.32256854 | 6.32256854 | < 1e-9 |
| $\beta_2$ | 19.04206847 | 19.04206847 | < 1e-9 |
| $\beta_3$ | -0.20682478 | -0.20682478 | < 1e-9 |
| SSE | 24037.651696 | 24037.651696 | < 1e-6 |
| RMSE | 15.504081 | 15.504081 | < 1e-6 |
| $R^2$ | 0.94830942 | 0.94830942 | < 1e-8 |
| $\operatorname{cond}_2(X)$ | 347.216931 | 347.216931 | < 1e-6 |
| $\sigma_i(X)$ | 572.522002, 80.371311, 32.809382, 1.648888 | idem | < 1e-6 |

**Correspondencia de funciones:**

| Operación | Octave | Eigen3 |
|---|---|---|
| Resolver $Ax=b$ | `A\b` | `A.partialPivLu().solve(b)` |
| LU con pivoteo | `[L,U,P]=lu(A)` | `PartialPivLU<MatrixXd>` |
| QR | `[Q,R]=qr(A)` | `HouseholderQR<MatrixXd>` |
| QR con pivoteo de columnas | `[Q,R,P]=qr(A)` | `ColPivHouseholderQr<MatrixXd>` |
| QR reducida (económica) | `[Q,R]=qr(A,0)` | `HouseholderQR` + `topRows(p)` |
| SVD | `[U,S,V]=svd(A)` | `JacobiSVD` / `BDCSVD` |
| Pseudoinversa | `pinv(A)` | `completeOrthogonalDecomposition().pseudoInverse()` |
| $\operatorname{cond}$ | `cond(A)` | `svd.singularValues()(0)/(...)(n-1)` |
| Cholesky | `chol(A)` | `LLT` / `LDLT` |

**Diferencias prácticas observadas:**

- Octave es **interpretado**: código mucho más corto, pero cada operación paga
  el costo del intérprete. Ideal para prototipar.
- Eigen es **header-only con templates y expresiones perezosas** (*expression
  templates*): `A*x - b` no crea temporales intermedios. Compilado con `-O2`
  aprovecha vectorización SIMD. Para matrices pequeñas (10×10, 100×4) Eigen es
  claramente más rápido; para matrices grandes ambos convergen porque los dos
  terminan llamando a BLAS/LAPACK.
- Eigen **no trae `cond()`**: hay que calcularlo desde la SVD.
- Los resultados coinciden hasta $\sim10^{-12}$ porque los algoritmos son los
  mismos; las últimas cifras difieren por el orden de las operaciones en punto
  flotante (que no es asociativo) y por el BLAS enlazado.

## Entregable 5 — Número de condición de X

| Cantidad | Valor |
|---|---|
| $\sigma(X)$ | 572.522002, 80.371311, 32.809382, **1.648888** |
| $\operatorname{cond}_2(X)$ | **347.216931** |
| $\operatorname{cond}_2(X^TX)$ | **120 559.597432** |
| $\operatorname{cond}_2(X)^2$ | 120 559.597432 **[OK]** |
| Dígitos perdidos con QR | 2.54 |
| Dígitos perdidos con ec. normales | 5.08 |
| $\operatorname{cond}_2(X)$ centrando columnas | **14.606127** |

Se verifica numéricamente la identidad $\operatorname{cond}(X^TX)=\operatorname{cond}(X)^2$.
**Formar las ecuaciones normales duplica los dígitos perdidos**: 2.5 → 5.1.

El origen del mal condicionamiento es la **escala dispar** de las columnas
(columna de unos vs. $T\sim50$), **no** colinealidad: centrar las columnas baja
$\operatorname{cond}$ de 347 a **14.6**, un factor 24, sin cambiar el modelo.
Esta es la medida preventiva más barata que existe en regresión.

## Diagnóstico adicional: el modelo físico

Un modelo lineal en $V$, $I$ y $T$ no es lo que dice la física. La potencia de
salida de un convertidor es $P=\eta\,V I$ — un **producto**, no una suma:

| | corr. con $P$ | RMSE | $R^2$ |
|---|---|---|---|
| Modelo lineal $\beta_0+\beta_1V+\beta_2I+\beta_3T$ | — | 15.5041 W | 0.94831 |
| **Modelo físico $a+b\,(VI)$** | **0.998434** | **3.8145 W** | **0.99687** |

Con **un solo regresor** ($V\!\cdot\!I$) el RMSE cae de 15.5 W a **3.8 W** (4×
mejor) y $R^2$ sube a 0.9969. El coeficiente $b=0.9105$ es directamente la
**eficiencia media del convertidor (≈91 %)**, y $a=-1.19$ W son las pérdidas en
vacío. El modelo lineal del enunciado es una **linealización** válida alrededor
del punto de operación medio, y es exactamente por eso que su RMSE residual
(11.5 % de la media) no es ruido sino curvatura no modelada.

## Efecto del número de muestras

| $n$ | $\beta_0$ | $\beta_1(V)$ | $\beta_2(I)$ | $\beta_3(T)$ | RMSE | $\operatorname{cond}(X)$ |
|---|---|---|---|---|---|---|
| 10 | -137.5630 | 6.7401 | 18.7333 | -0.0830 | 13.8115 | 236.59 |
| 20 | -136.0298 | 6.4210 | 19.0120 | -0.0454 | 15.1012 | 287.85 |
| 50 | -123.3646 | 6.2820 | 19.0565 | -0.2228 | 15.0082 | 336.00 |
| 100 | -123.4183 | 6.3226 | 19.0421 | -0.2068 | 15.5041 | 347.22 |

Errores estándar de $\beta_1$: 0.944 ($n$=10) → 0.637 (20) → 0.390 (50) → 0.268 (100).
La reducción sigue aproximadamente $1/\sqrt{n}$: multiplicar $n$ por 10 reduce el
error estándar en un factor 3.5, cercano al $\sqrt{10}=3.16$ teórico.

---
# PARTE 3 — Preguntas de análisis

## 1. ¿Qué variable tiene mayor impacto en la potencia?

**La corriente de carga $I$**, de forma inequívoca y por tres criterios
independientes que coinciden:

| Criterio | $V$ | $I$ | $T$ |
|---|---|---|---|
| $\beta$ estandarizado | 0.549 | **0.814** | -0.044 |
| Influencia relativa | 39.0 % | **57.9 %** | 3.1 % |
| Efecto en todo el rango | +124.1 W | **+186.4 W** | -10.2 W |
| Estadístico $t$ | 23.6 | **34.7** | -1.9 |
| Correlación con $P$ | 0.523 | **0.802** | -0.184 |

$\beta_2=19.04$ es mayor que $\beta_1=6.32$, pero eso no significa nada por sí
solo: son unidades distintas (W/A vs. W/V). Si la corriente se midiera en
miliamperios, $\beta_2$ valdría 0.019 y el modelo sería idéntico. Solo tras
estandarizar por la dispersión de cada variable la comparación es válida.

Físicamente tiene todo el sentido: $P=VI$, y en estos datos la corriente varía
proporcionalmente más que el voltaje ($s_I/\bar I = 42\%$ frente a
$s_V/\bar V = 28\%$).

## 2. ¿Cómo influye la temperatura en el rendimiento?

Negativamente, pero de forma estadísticamente no significativa en estos datos.

- $\beta_3=-0.2068$ W/°C: por cada grado de aumento del disipador, la potencia
  de salida baja ~0.21 W.
- Sobre todo el rango medido (25 → 74.5 °C): **-10.2 W**, apenas el 3.1 % de la
  influencia total.
- $t=-1.886$, $p=0.062 > 0.05$: **no se rechaza $H_0:\beta_3=0$** al 5 % de
  significancia. El intervalo de confianza al 95 % incluye el cero.

El signo es físicamente correcto y esperable: al subir la temperatura
aumenta la resistencia $R_{DS(on)}$ de los MOSFET (coeficiente térmico positivo,
típicamente +0.4 %/°C en silicio), crecen las pérdidas por conducción
$P_{\text{cond}}=I^2R_{DS(on)}$ y baja la eficiencia. También aumentan las
pérdidas en el cobre del inductor.

La magnitud es pequeña y el efecto queda enmascarado por
el ruido. Para caracterizar el efecto térmico adecuadamente haría falta (i) más
muestras, (ii) barrer $T$ manteniendo $V$ e $I$ fijos, o (iii) modelar la
interacción real $P \sim \eta(T)\,VI$, ya que la temperatura afecta la
eficiencia (un factor multiplicativo), no la potencia de forma aditiva. Un
término de interacción $T\cdot V I$ sería más apropiado que el término aditivo
$\beta_3 T$.

## 3. ¿Qué ocurre al aumentar el número de muestras?

Cuatro efectos distintos, que conviene no confundir:

**(a) Los coeficientes se estabilizan y el sesgo desaparece.** De $n$=10 a
$n$=100, $\beta_0$ pasa de -137.6 a -123.4 y $\beta_3$ de -0.083 a -0.207. Con
10 muestras y 4 parámetros solo quedan 6 grados de libertad: la estimación está
dominada por el azar del muestreo.

**(b) Los errores estándar decrecen como $1/\sqrt{n}$.** Para $\beta_1$:
0.944 → 0.637 → 0.390 → 0.268. **Cuadruplicar $n$ solo reduce el error a la
mitad**: es un retorno decreciente, y por eso pasar de 100 a 1000 muestras
compra menos de lo que la gente espera.

**(c) El RMSE NO baja; converge.** 13.81 → 15.10 → 15.01 → 15.50 W. Esto
sorprende a mucha gente pero es correcto: el RMSE estima el **ruido real del
proceso**, que no depende de $n$. Con pocas muestras el RMSE aparece
*artificialmente bajo* porque el modelo sobreajusta (con $n=p$ el residual sería
exactamente cero y el modelo no valdría nada). Lo que mejora con $n$ no es el
ajuste a los datos vistos, sino la capacidad predictiva sobre datos nuevos.

**(d) El condicionamiento se estabiliza.** $\operatorname{cond}(X)$: 236.6 →
287.8 → 336.0 → 347.2, tendiendo a un valor asintótico determinado por el
diseño experimental (los rangos y las escalas de las variables), no por $n$.
Más muestras no arreglan un problema mal condicionado; solo un mejor diseño
o un escalado adecuado lo hacen.

**(e) Costo computacional.** QR es $O(np^2)$: **lineal** en $n$. Duplicar las
muestras duplica el tiempo. Con $p$ fijo y pequeño esto escala sin problema a
millones de filas (y admite formulación incremental / por bloques).

## 4. ¿Por qué QR es preferible a las ecuaciones normales?

$$\operatorname{cond}_2(X^TX)=\operatorname{cond}_2(X)^2$$

En nuestros datos: $347.216931^2 = 120\,559.597432$, exactamente el valor
calculado. **Formar $X^TX$ duplica los dígitos perdidos** (2.54 a 5.08).

Las razones, ordenadas por importancia:

1. **Amplificación del condicionamiento.** El error relativo de las ecuaciones
   normales escala con $\operatorname{cond}(X)^2\varepsilon$; el de QR, con
   $\operatorname{cond}(X)\varepsilon$. Si $\operatorname{cond}(X)=10^8$ —
   perfectamente posible en un ajuste polinómico o con variables de escalas muy
   distintas — las ecuaciones normales dan $10^{16}\cdot10^{-16}=1$: **cero
   cifras correctas**, mientras QR conserva 8.

2. **Pérdida de información al formar el producto.** El ejemplo canónico
   (Läuchli): con $\epsilon=10^{-9}$,
   $X=\begin{bmatrix}1&1\\ \epsilon&0\\0&\epsilon\end{bmatrix}$ tiene rango 2,
   pero $X^TX=\begin{bmatrix}1+\epsilon^2&1\\1&1+\epsilon^2\end{bmatrix}$ se
   redondea a $\begin{bmatrix}1&1\\1&1\end{bmatrix}$ porque
   $1+\epsilon^2 = 1$ en doble precisión: **la matriz se vuelve singular y la
   información del rango se destruye irreversiblemente.**

3. **Invariancia ortogonal.** $\|Q^Tv\|_2=\|v\|_2$: las transformaciones de
   Householder no amplifican el error ($\operatorname{cond}_2(Q)=1$). En
   este taller se verificó $\|Q^TQ-I\|_F\approx10^{-15}$.

4. **QR no requiere rango completo** si se usa pivoteo de columnas
   (`ColPivHouseholderQr`), que además detecta la deficiencia de rango.
   Cholesky sobre $X^TX$ simplemente falla.

5. **Estabilidad demostrada.** QR con Householder es *backward stable*: el
   $\hat\beta$ calculado es la solución exacta de un problema con $X+\delta X$,
   donde $\|\delta X\|/\|X\|=O(\varepsilon)$ — la mejor garantía posible.

**¿Cuándo sí usar ecuaciones normales?** Cuando $p$ es muy pequeño, $n\gg p$,
$X$ está bien condicionada y el costo importa: son ~2× más baratas
($np^2$ vs. $2np^2$) y $X^TX$ se puede acumular por bloques o en streaming. En
este taller, con $\operatorname{cond}(X)=347$, dan el mismo resultado —
pero con un error 2–3 órdenes de magnitud mayor, que es el aviso de lo que pasa
cuando el condicionamiento crece.

**Y la SVD, ¿cuándo?** Cuando el rango es dudoso, hay colinealidad severa, o se
necesita regularizar. Cuesta ~4× lo que QR, pero es la única que diagnostica.

## 5. ¿Cómo se interpreta físicamente cada coeficiente?

$$\hat P = -123.4183 + 6.3226\,V + 19.0421\,I - 0.2068\,T$$

**$\beta_0 = -123.42$ W (intercepto).**
Matemáticamente: potencia predicha con $V=I=T=0$. **No tiene interpretación
física directa** — es una extrapolación muy fuera del rango medido
($V\in[12.1, 31.7]$, $I\in[2.07, 11.86]$), y una potencia negativa es
físicamente imposible. Es un **parámetro de ajuste** que absorbe la curvatura de
la superficie real $P=\eta VI$ al aproximarla por un plano. Su valor grande y
negativo es, de hecho, *síntoma* de que el modelo lineal está forzando un plano
sobre una superficie bilineal. Si se centran las variables, el intercepto pasa a
ser $\bar P = 134.27$ W y sí adquiere sentido: la potencia en el punto de
operación medio.

**$\beta_1 = +6.3226$ W/V (voltaje).**
Derivada parcial $\partial P/\partial V$ a corriente y temperatura constantes.
Del modelo físico $P=\eta VI$ se tiene $\partial P/\partial V=\eta I$, que
evaluado en la corriente media da $0.91 \times 6.98 = 6.35$ W/V — **prácticamente
idéntico al 6.32 estimado**. El coeficiente es por tanto la *eficiencia por la
corriente media*, y su validez está limitada al entorno de $\bar I$.

**$\beta_2 = +19.0421$ W/A (corriente).**
$\partial P/\partial I = \eta V \approx 0.91\times 21.40 = 19.49$ W/A frente a
los **19.04** estimados. La coincidencia vuelve a confirmar la lectura: es la
*eficiencia por el voltaje medio*. Es el coeficiente dominante porque el voltaje
medio (21.4 V) es numéricamente mayor que la corriente media (6.98 A).

**$\beta_3 = -0.2068$ W/°C (temperatura).**
Degradación térmica: pérdida de ~0.21 W de salida por cada grado del disipador.
Recoge el aumento de $R_{DS(on)}$ de los interruptores y de la resistencia del
cobre con la temperatura. Es el único coeficiente no significativo
($p=0.062$), así que debe reportarse como *tendencia observada consistente con
la física*, no como efecto establecido. Físicamente su lugar correcto sería
multiplicativo, $\eta = \eta_0(1-\alpha T)$, no aditivo.

**Verificación cruzada de conjunto.** El modelo alternativo $P = a + b(VI)$ da
$b = 0.9105$ con $R^2=0.9969$. Que $\beta_1 \approx b\,\bar I$ y
$\beta_2 \approx b\,\bar V$ se cumplan simultáneamente con 2 cifras de acuerdo no
es casualidad: es la confirmación de que los coeficientes lineales son las
**derivadas parciales de la superficie física evaluadas en el punto de operación
medio**, exactamente lo que predice una expansión de Taylor de primer orden:

$$P \approx \eta\bar V\bar I + \eta\bar I(V-\bar V) + \eta\bar V(I-\bar I) + \mathcal{O}(\Delta^2)$$

Los coeficientes de una regresión lineal multivariable **no son constantes
universales del sistema**: son válidos solo en el entorno donde se tomaron los
datos. Ese es el mensaje de ingeniería de todo el problema.

---

## Referencias

1. Golub, G. H. y Van Loan, C. F. (2013). *Matrix Computations*, 4ª ed. Johns Hopkins University Press.
2. Trefethen, L. N. y Bau, D. (1997). *Numerical Linear Algebra*. SIAM.
3. Higham, N. J. (2002). *Accuracy and Stability of Numerical Algorithms*, 2ª ed. SIAM.
4. Guennebaud, G., Jacob, B. et al. (2010). *Eigen v3*. https://eigen.tuxfamily.org
5. Eaton, J. W. et al. *GNU Octave Reference Manual*. https://octave.org/doc/
6. Erickson, R. W. y Maksimović, D. (2020). *Fundamentals of Power Electronics*, 3ª ed. Springer.

