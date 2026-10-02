# Taller: Sistema masa–resorte–amortiguador (EDO de segundo orden)

**Métodos Computacionales en Ingeniería (MCEI_M)**
Escuela Colombiana de Ingeniería Julio Garavito
Autor: Carlos Castillo

## Problema

El movimiento de un sistema acoplado masa–resorte está descrito por la EDO

$$m\,\ddot x + c\,\dot x + k\,x = 0$$

con $m = 20$ kg, $k = 20$ N/m, $x(0) = 1$ m y $\dot x(0) = 0$. Se resuelve numéricamente en $0 < t < 15$ s para tres coeficientes de amortiguamiento y se grafican los tres desplazamientos sobre la misma figura.

| $c$ (N·s/m) | $\zeta = c / 2\sqrt{mk}$ | Régimen |
|---|---|---|
| 5 | 0.125 | Subamortiguado |
| 40 | 1 | Crítico |
| 200 | 5 | Sobreamortiguado |

## Método

1. **Reducción a primer orden.** Con $\mathbf{y} = [x,\ \dot x]^T$:
   $\dot y_1 = y_2,\quad \dot y_2 = -\frac{c}{m}y_2 - \frac{k}{m}y_1$.
2. **Runge–Kutta clásico de 4.º orden (RK4)**, implementado a mano con paso fijo $h = 0.01$ s.
3. **Validación** contra la solución analítica de cada régimen y contra `scipy.integrate.solve_ivp` (RK45 adaptativo).
4. **Verificación del orden de convergencia** de RK4 variando $h$.

## Archivos

| Archivo | Descripción |
|---|---|
| `masa_resorte_CarlosCastillo.ipynb` | Cuadernillo con planteamiento, implementación de RK4, validación, gráficas y conclusiones |
| `masa_resorte_desplazamiento.png` | Desplazamiento vs. tiempo para los tres coeficientes de amortiguamiento |

## Ejecución

```bash
pip install numpy scipy matplotlib jupyter
jupyter lab masa_resorte_CarlosCastillo.ipynb
```

También se puede abrir directamente en VS Code con la extensión de Jupyter. La figura se guarda automáticamente como `masa_resorte_desplazamiento.png`.

## Resultados

![Desplazamiento vs tiempo](masa_resorte_desplazamiento.png)

| $c$ (N·s/m) | Error máx. RK4 ($h=0.01$ s) | Error máx. `solve_ivp` (RK45) |
|---|---|---|
| 5 | 2.5e-10 m | 4.3e-9 m |
| 40 | 9.7e-11 m | 6.3e-10 m |
| 200 | 3.3e-9 m | 1.2e-10 m |

**Convergencia de RK4** (caso $c = 5$): al dividir $h$ por 2, el error se divide por ≈ 16, lo que confirma el orden $p \approx 4$.

| $h$ (s) | Error máx. | Orden observado |
|---|---|---|
| 0.4 | 6.2e-4 | – |
| 0.2 | 4.0e-5 | 3.94 |
| 0.1 | 2.5e-6 | 4.02 |
| 0.05 | 1.5e-7 | 4.01 |
| 0.025 | 9.7e-9 | 4.00 |

## Conclusiones

- **Subamortiguado ($c=5$):** oscila con periodo $T_d \approx 6.33$ s y amplitud que decae como $e^{-0.125t}$. A los 15 s todavía no se ha detenido.
- **Crítico ($c=40$):** vuelve al equilibrio en el menor tiempo posible sin oscilar.
- **Sobreamortiguado ($c=200$):** no oscila, pero retorna mucho más lento que el crítico ($x(15) \approx 0.22$ m), porque domina el modo lento $e^{-0.10t}$. Más amortiguamiento no implica un retorno más rápido.
- RK4 con $h = 0.01$ s reproduce la solución exacta con errores menores a $10^{-8}$ m. En el caso sobreamortiguado, las dos escalas de tiempo ($|r_2/r_1| \approx 98$) introducen una rigidez leve: la raíz rápida limita el paso estable de un método explícito ($h \lesssim 0.28$ s para RK4).
