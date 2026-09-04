---------------------------------------------
Resultados 
---------------------------------------------
**Preguntas de Exploración:**
1. ¿Cómo identificar visualmente la existencia de una raíz?
R/: Una raíz de una función f(x)es un valor para el cual: f(x)=0. Visualmente, una raíz puede identificarse observando el punto donde la gráfica de la función intersecta o toca el eje x. Otra forma es observar si la función cambia de signo dentro de un intervalo. Si: f(a)⋅f(b)<0 y la función es continua en [aⓜ,b], entonces existe al menos una raíz dentro de ese intervalo, de acuerdo con el Teorema del Valor Intermedio.

2. ¿Qué significa que una raíz esté acotada?
R/: Significa que podemos establecer un intervalo [aⓜ,b]dentro del cual se encuentra la raíz. Por ejemplo, si f(a)<0 "y"  f(b)>0 se puede decir que bajo las condiciones de continuidad, que existe al menos una raíz entre ay b. Este concepto es fundamental en métodos cerrados, como el método de bisección, porque estos métodos mantienen la raíz dentro de un intervalo durante el proceso iterativo.

3. ¿Por qué algunos métodos requieren derivadas?
R/: Algunos métodos utilizan la pendiente de la función para determinar la dirección y magnitud del siguiente desplazamiento hacia la raíz. Por ejemplo en el método de Newton-Raphson se utiliza la derivada f^' (x) para conocer la pendiente de la función en el punto actual y aproximar la recta tangente donde se intersectará el eje x.
La ventaja es que puede producir una convergencia muy rápida cerca de la raíz. Sin embargo, requiere calcular la derivada y puede presentar problemas cuando esta es cero o cercana a cero.

4. ¿Qué ventajas tiene un método abierto frente a uno cerrado?
R/: Los métodos abiertos, como Newton-Raphson y el método de la secante, no necesitan mantener la raíz dentro de un intervalo previamente acotado. Su principal ventaja es que pueden presentar una convergencia mucho más rápida que métodos cerrados como bisección. Sin embargo, esta ventaja no siempre garantiza la convergencia. Por otro lado, los métodos cerrados son generalmente más robustos porque mantienen la raíz dentro de un intervalo que la contiene.

5. ¿Existe garantía de convergencia en todos los casos?
R/: No. Por ejemplo:
	Bisección: tiene una garantía de convergencia si la función es continua y existe un cambio de signo en el intervalo. 
	Newton-Raphson: puede no converger si la aproximación inicial es inadecuada, si f^' (x)es cero o si la función presenta determinadas características. 
	Secante: tampoco garantiza convergencia y puede divergir dependiendo de los valores iniciales. 
	Regula Falsi: generalmente es robusto porque conserva el intervalo acotado, aunque puede presentar convergencia lenta en algunos casos. 
Por esta razón, la elección del método depende de las características de la función y del problema.

**Resultados Ejercicio función f(x) = x^3 - 5x + 1**
Se evaluacron la funcion en los intervalos [-3,2], encontrando 3 raices:
/n
**Primera Raiz**
Método    | raíz | iteraciones | Observaciones |
Bisección | -2.3300587460 | 26 |  converge|
Falsa Posición | -2.3300587396 | 8 | converge |
Brent | -2.3300587396 | 6 | converge |
Newton | -2.3300587396 | 4 | converge |
Secante | -2.3300587396 | 5 | converge |
Steffenson | -2.3300587396 | 5 | converge |
/n 
**Segunda Raiz**
Método    | raíz | iteraciones | Observaciones |
Bisección | 0.2016396755 | 29 | converge |
Falsa Posición | 0.2016396757 | 8 | converge |
Brent | 0.2016396757 | 6 | converge |
Newton | 0.2016396757 | 4 | converge |
Secante | 0.2016396757 | 5 | converge |
Steffenson | 0.2016396757 | 4 | converge |
/n
**Tercera Raiz**
Método    | raíz | iteraciones | Observaciones |
Bisección | 2.1284190640 | 26 | converge |
Falsa Posición | 2.1284190638 | 25 | converge |
Brent | 2.1284190638 | 6 | converge |
Newton | 2.1284190638 | 5 | converge |
Secante | 2.1284190638 | 6 | converge |
Steffenson | 2.1284190638 | 6 | converge |

**Preguntas de Discusión**

1. ¿Qué método converge en menos iteraciones?
R/: El método de Newton presenta la menor cantidad de iteraciones en los experimentos realizados. Para las tres raíces obtuvo:

- Raíz 1: 4 iteraciones
- Raíz 2: 4 iteraciones
 -Raíz 3: 5 iteraciones

Por tanto, en este experimento fue el método más rápido en términos de iteraciones. La secante y Steffensen requirieron generalmente 5–6 iteraciones, mientras que Brent necesitó 6 iteraciones

2. ¿Qué método presenta mayor robustez?
R/: El método que presenta mayor robustez es la bisección, debido a que mantiene la raíz dentro de un intervalo y reduce progresivamente su amplitud. En los tres casos siempre convergió, aunque necesitó más iteraciones: 26, 29 y 26, respectivamente. Sin embargo, Brent ofrece una excelente relación entre robustez y rapidez, ya que también trabajó con intervalos acotados y convergió en solamente 6 iteraciones en los tres experimentos.

3. ¿Qué ocurre cuando el valor inicial está lejos de la raíz?
R/: Cuando se utiliza un método abierto, como Newton, Secante o Steffensen, un valor inicial alejado de la raíz puede provocar que el método necesite más iteraciones, presente oscilaciones o incluso no converja, dependiendo de la función.
En los resultados analizados, los valores iniciales utilizados fueron relativamente cercanos a las raíces, por lo que no se observa un caso de divergencia. Por ejemplo, Newton partiendo de \(x=2.5\) alcanzó la raíz \(2.1284190638\) en 5 iteraciones.
Por otro lado, los métodos cerrados como bisección, falsa posición y Brent trabajan con intervalos que contienen la raíz, por lo que son menos sensibles a la elección de un punto inicial individual.

4. ¿Cuál método recomendaría para problemas de ingeniería?
R/: Para problemas de ingeniería recomendaría Brent cuando se disponga de un intervalo que contenga la raíz, ya que convergio el resultado en pocas iteraciones y sin necesidad de la derivada de la función. Y newton en los demas casos.

5. ¿Existe una relación entre costo computacional y velocidad de convergencia?
R/: Sí. Existe un compromiso entre costo computacional por iteración y cantidad de iteraciones. 
Por ejemplo, Newton converge rápidamente, pero en cada iteración requiere evaluar tanto la función como su derivada, en este caso necesitó solamente 4–5 iteraciones.

Por otro lado, La bisección, tiene operaciones sencillas y bajo costo por iteración, pero requiere muchas más iteraciones: hasta 29 en los experimentos realizados.

Por ello, el método con menos iteraciones no necesariamente tiene siempre el menor costo total. Para seleccionar un método de ingeniería se debe considerar no solo el número de iteraciones, sino también el costo de evaluar la función, calcular derivadas, memoria, estabilidad y garantía de convergencia.
---------------------------------------------
