# Guía para Dummies: Métodos Cuasi-Newtonianos, BFGS y Newton-CG

> **Materia:** Métodos Numéricos II (`DAT-252`) — Carrera de Informática, UMSA  
> **Propósito de esta guía:** Entender **de verdad** qué hacen BFGS y Newton-CG sin ahogarse en derivadas parciales, y saber defender cada fórmula ante cualquier pregunta del docente.

---

## 1. La Analogía Maestra: Bajando una Montaña con Neblina

Imagina que estás perdido en una cordillera de noche y con neblina espesa. Tu meta es llegar al **punto más bajo del valle** (eso es **minimizar una función** $f(x)$). Como no ves el mapa completo, solo puedes medir cosas en el punto $x_k$ donde estás parado ahora mismo.

Hay **4 tipos de excursionistas** (los 4 grandes métodos de optimización):

### 1. Descenso de Gradiente: *El excursionista miope*
- **Qué hace:** Siente con los pies hacia dónde baja el suelo ($\nabla f$) y da un paso en esa dirección.
- **Problema:** Si estás en un cañón muy largo y angosto (como un tobogán curvo: la *función de Rosenbrock*), la pendiente más empinada siempre apunta hacia la pared de enfrente, no a lo largo del cañón. Resultado: rebota de pared en pared en **zig-zag** durante miles de pasos.

### 2. Método de Newton Exacto: *El millonario con dron escáner láser 3D*
- **Qué hace:** En cada paso saca un dron que escanea la **curvatura completa** del terreno en todas las direcciones (la **Matriz Hessiana** $\nabla^2 f$). Con eso dibuja un «cuenco parabólico» perfecto y salta directo al fondo en muy pocos pasos (convergencia cuadrática).
- **Dos problemas enormes:**
  1. **Carísimo:** Si el terreno tiene $n$ variables, medir la curvatura requiere $n^2$ derivadas y resolver el salto cuesta $O(n^3)$ operaciones. Con $n = 10{,}000$, solo guardar la matriz ocupa ~800 MB y factorizarla toma minutos por paso.
  2. **Peligroso en zonas raras (no convexas):** Si estás cerca de una cima o un **paso de montaña (punto silla)** donde el terreno se curva hacia abajo en una dirección, el escáner confunde «fondo» con «cima» y te lanza **cuesta arriba** o directo hacia la trampa del punto silla.

### 3. Métodos Cuasi-Newtonianos y BFGS: *El excursionista con memoria y cuaderno*
- **Qué hace:** ¡No compra ningún dron! Simplemente da un paso ($s_k$) y anota en su cuaderno **cuánto cambió la pendiente del suelo** entre donde estaba y donde llegó ($y_k$).
  - Si caminó mucho y la pendiente casi no cambió $\rightarrow$ el suelo es plano (puedo dar pasos grandes).
  - Si caminó poquito y la pendiente cambió bruscamente $\rightarrow$ el suelo es muy curvo (pasos cortos y precisos).
- **La magia de BFGS:** En cada paso actualiza su «mapa aproximado» de la curvatura inversa ($H_k$) usando solo restas y multiplicaciones de vectores ($O(n^2)$, nunca invierte matrices ni calcula segundas derivadas). Además, se asegura de que su mapa **siempre sea un cuenco hacia arriba** (definido positivo), así que **nunca te manda cuesta arriba**.

### 4. Método de Newton-CG: *El excursionista en 1 millón de dimensiones*
- **Qué problema resuelve:** Cuando $n = 1{,}000{,}000$ (redes neuronales, física computacional, imágenes), **ni siquiera cabe el cuaderno de BFGS en la RAM** ($n \times n$ ocuparía 8 Terabytes).
- **Qué hace:** En vez de dibujar el mapa completo de la curvatura, le hace **unas pocas preguntas rápidas** al terreno en direcciones inteligentes (**Gradiente Conjugado — CG**): *"si me muevo en dirección $d$, ¿cómo cambia la pendiente?"* (producto Hessiana-vector $\nabla^2 f \cdot d$).
  - **Sin matriz (*Hessian-Free*):** Esa pregunta se responde con una resta de dos gradientes, usando solo $O(n)$ de RAM (4 vectores).
  - **Sin perfeccionismo inútil (*Newton Inexacto*):** Cuando estás lejos del fondo, no pierdes tiempo calculando el salto con 16 decimales; paras CG temprano según una tolerancia $\eta_k$.
  - **Freno de emergencia (curvatura negativa):** Si en una pregunta detecta que el terreno se curva al revés ($d^T \nabla^2 f\, d \le 0$), **corta CG de inmediato** y usa la mejor dirección segura que ya tenía.

---

## 2. Diccionario «Sin Miedo» (Traduciendo los Símbolos)

| Símbolo | Nombre formal | Traducción a español humano |
| :---: | :--- | :--- |
| $f(x)$ | Función objetivo | La **altura** en la montaña que quieres minimizar. |
| $\nabla f(x_k)$ | Gradiente | Flecha que apunta hacia la **subida más empinada**. Por eso $-\nabla f$ apunta cuesta abajo. |
| $\nabla^2 f(x_k)$ | Matriz Hessiana | La **curvatura exacta** del terreno en todas las direcciones (matriz $n \times n$ simétrica). |
| $p_k$ | Dirección de búsqueda | La **brújula**: hacia dónde decides caminar en el paso $k$. |
| $\alpha_k$ | Tamaño de paso (*step length*) | **Qué tan largo** das el paso en la dirección $p_k$ ($x_{k+1} = x_k + \alpha_k p_k$). |
| $s_k$ | Paso efectivo | $s_k = x_{k+1} - x_k$: cuánto te moviste realmente. |
| $y_k$ | Cambio de gradiente | $y_k = \nabla f(x_{k+1}) - \nabla f(x_k)$: cuánto cambió la pendiente tras dar ese paso. |
| $B_k$ | Aproximación de la Hessiana | Tu **estimación** de la curvatura $\nabla^2 f(x_k)$. |
| $H_k$ | Inversa de $B_k$ ($H_k = B_k^{-1}$) | El «mapa listo para multiplicar»: te da el salto directo $p_k = -H_k \nabla f_k$ sin resolver sistemas lineales. |
| $\rho_k$ | Escalar $\frac{1}{y_k^T s_k}$ | Un factor de normalización que aparece en BFGS. Para que no explote, necesitamos $s_k^T y_k > 0$. |
| $A \succ 0$ | Matriz definida positiva (SPD) | Significa que la superficie tiene forma de **tazón hacia arriba** (tiene fondo único). Garantiza que $-H_k \nabla f_k$ siempre baja. |
| $\eta_k$ | Secuencia de forzamiento | Qué tan «exigente» es Newton-CG en la iteración $k$. Grande ($\approx 0.5$) lejos del óptimo; diminuto ($\to 0$) cerca del óptimo. |
| **Lineal / Superlineal / Cuadrática** | Velocidad de convergencia | **Lineal:** ganas una fracción fija de dígitos por paso (lento). **Superlineal (BFGS):** cada paso es proporcionalmente mejor que el anterior. **Cuadrática (Newton):** el número de decimales correctos se **duplica** en cada paso ($1 \to 2 \to 4 \to 8 \to 16$). |

---

## 3. Entendiendo BFGS en 3 Ideas (Sin Memorizar la Fórmula a Ciegas)

### Idea 1: La Ecuación Secante (1D $\rightarrow$ $n$D)
En cálculo de 1 variable, si no quieres calcular la segunda derivada $f''(x)$, la aproximas dividiendo el cambio de la primera derivada entre el cambio de $x$:
$$f''(x_{k+1}) \approx \frac{f'(x_{k+1}) - f'(x_k)}{x_{k+1} - x_k} = \frac{y_k}{s_k}$$
En $n$ dimensiones **no se pueden dividir vectores**, así que pasamos $s_k$ multiplicando al otro lado:
$$B_{k+1} s_k = y_k \quad \Longleftrightarrow \quad H_{k+1} y_k = s_k$$
¡Eso es toda la **ecuación secante**! Le pide al nuevo mapa $B_{k+1}$ que sea coherente con el último paso que acabas de dar.

### Idea 2: ¿Por qué necesitamos $s_k^T y_k > 0$? (Condición de Curvatura)
Multiplica la ecuación secante por $s_k^T$:
$$s_k^T B_{k+1} s_k = s_k^T y_k$$
Para que $B_{k+1}$ represente un «tazón hacia arriba» (definida positiva), el lado izquierdo **tiene** que ser positivo. Por lo tanto, **es obligatorio que $s_k^T y_k > 0$**.
- **¿Cómo lo garantizamos siempre?** Eligiendo el largo de paso $\alpha_k$ con las **Condiciones de Wolfe**:
  1. **Armijo:** «Asegúrate de bajar lo suficiente» (evita pasos gigantes que empeoran la altura).
  2. **Curvatura:** «Asegúrate de que la pendiente se suavizó un poco» (evita pasos diminutos y fuerza automáticamente $s_k^T y_k > 0$).

### Idea 3: La Fórmula de BFGS (y por qué actualiza $H_k$ en vez de $B_k$)
Hay infinitas matrices que cumplen $H_{k+1} y_k = s_k$. BFGS elige la que **menos cambia** tu mapa anterior $H_k$ (para no olvidar lo aprendido en pasos previos) y actualiza directamente la inversa $H_k$ usando el truco de **Sherman–Morrison–Woodbury**:

$$H_{k+1} = (I - \rho_k s_k y_k^T)\, H_k\, (I - \rho_k y_k s_k^T) + \rho_k s_k s_k^T, \qquad \rho_k = \frac{1}{y_k^T s_k}$$

**Traducción humana de esa fórmula:**
1. $(I - \rho_k s_k y_k^T) H_k (I - \rho_k y_k s_k^T)$: «Toma tu mapa viejo $H_k$ y **limpia** la parte que contradecía lo que acabas de medir en la dirección $y_k$».
2. $+ \rho_k s_k s_k^T$: «**Pega** en ese lugar la curvatura fresca que acabas de observar ($s_k$ con $y_k$)».
3. **Gran premio:** Como actualizamos $H_k$ directamente, calcular $p_k = -H_k \nabla f_k$ es una simple multiplicación matriz-vector en **$O(n^2)$**, ¡nunca resolvemos un sistema lineal en $O(n^3)$!

---

## 4. Mini-Ejemplo de BFGS Hecho a Mano (1 Iteración Completa)

Supongamos que queremos minimizar una elipse sencilla:
$$f(x_1, x_2) = x_1^2 + 5 x_2^2 \qquad \Longrightarrow \qquad \nabla f(x) = \begin{bmatrix} 2x_1 \\ 10x_2 \end{bmatrix}, \quad \nabla^2 f_{\text{real}} = \begin{bmatrix} 2 & 0 \\ 0 & 10 \end{bmatrix}, \quad (\nabla^2 f)^{-1}_{\text{real}} = \begin{bmatrix} 0.5 & 0 \\ 0 & 0.1 \end{bmatrix}$$

Arrancamos en $x_0 = \begin{bmatrix} 1 \\ 1 \end{bmatrix}$ con $H_0 = I = \begin{bmatrix} 1 & 0 \\ 0 & 1 \end{bmatrix}$ (no sabemos nada de la curvatura todavía):

1. **Gradiente inicial:** $\nabla f(x_0) = \begin{bmatrix} 2 \\ 10 \end{bmatrix}$.
2. **Dirección inicial:** $p_0 = -H_0 \nabla f(x_0) = \begin{bmatrix} -2 \\ -10 \end{bmatrix}$ (el primer paso de BFGS siempre coincide con Descenso de Gradiente).
3. **Búsqueda de línea ($\alpha_0$):** El paso exacto a lo largo de $p_0$ da $\alpha_0 = \frac{2^2 + 10^2}{2(2^2) + 10(10^2)} = \frac{104}{1008} \approx 0.1032$.
4. **Nuevo punto:** $x_1 = x_0 + \alpha_0 p_0 \approx \begin{bmatrix} 0.7937 \\ -0.0317 \end{bmatrix}$.
5. **Medimos $s_0$ y $y_0$:**
   $$s_0 = x_1 - x_0 \approx \begin{bmatrix} -0.2063 \\ -1.0317 \end{bmatrix}, \qquad \nabla f(x_1) \approx \begin{bmatrix} 1.5873 \\ -0.3175 \end{bmatrix}, \qquad y_0 = \nabla f(x_1) - \nabla f(x_0) \approx \begin{bmatrix} -0.4127 \\ -10.3175 \end{bmatrix}$$
6. **Verificamos curvatura:** $s_0^T y_0 \approx (-0.2063)(-0.4127) + (-1.0317)(-10.3175) \approx 10.73 > 0$. ¡Perfecto! $\rho_0 = 1/10.73 \approx 0.0932$.
7. **Actualizamos $H_1$ con la fórmula BFGS:**
   $$H_1 \approx \begin{bmatrix} 0.984 & -0.038 \\ -0.038 & 0.102 \end{bmatrix}$$
   **¡Mira el número en la esquina inferior derecha ($0.102$)!** En un **solo paso**, sin calcular una sola segunda derivada, BFGS ya descubrió que en $x_2$ la inversa de la curvatura es $\approx 0.1$ (el valor real exacto es $0.1$). En el siguiente paso ($k=2$), BFGS acierta directo al mínimo $(0, 0)$!

---

## 5. Entendiendo Newton-CG en 3 Ideas

Si ya expusieron **Newton-Krylov** en su primera exposición, ¡ya entienden el 80 % de Newton-CG!
- En Newton-Krylov resolvían $F(x) = 0$ con Jacobiano $J(x)$ **no simétrico** $\implies$ usaban **GMRES**.
- En optimización resolvemos $\nabla f(x) = 0$, cuyo Jacobiano es la **Hessiana $\nabla^2 f(x)$**, que **SIEMPRE es simétrica** $\implies$ cambiamos GMRES por **Gradiente Conjugado (CG)**, que es mucho más barato en memoria (solo 4 vectores, $O(n)$).

Las 3 reglas de oro de Newton-CG:
1. **No seas perfeccionista lejos de la meta (Secuencia de forzamiento $\eta_k$):**
   No resuelvas $\nabla^2 f(x_k) p_k = -\nabla f_k$ con 16 decimales cuando estás lejos de $x^*$. Para CG en cuanto:
   $$\|\nabla^2 f(x_k) p_k + \nabla f(x_k)\| \le \eta_k \|\nabla f(x_k)\|$$
   Con $\eta_k = \min(0.5, \sqrt{\|\nabla f_k\|})$ obtienes convergencia **superlineal**; con $\eta_k = \min(0.5, \|\nabla f_k\|)$ recuperas convergencia **cuadrática**.
2. **No construyas la matriz Hessiana (*Hessian-Free*):**
   CG solo pide el vector $\nabla^2 f(x_k) d$. Si no quieres programar la Hessiana, lo aproximas con una resta de gradientes:
   $$\nabla^2 f(x_k)\, d \approx \frac{\nabla f(x_k + h d) - \nabla f(x_k)}{h}, \qquad h \approx \frac{1.5 \times 10^{-8}}{\|d\|}$$
3. **Freno automático ante curvatura negativa:**
   En cada paso $j$ dentro de CG se calcula el número $\kappa_j = d_j^T \nabla^2 f(x_k) d_j$.
   - Si $\kappa_j > 0$: la curvatura en esa dirección es un tazón hacia arriba. CG sigue normal.
   - Si $\kappa_j \le 0$: ¡ALERTA! Pisaste una zona cóncava o un punto silla. Si estás en el primer paso ($j=0$), sales usando descenso de gradiente $p_k = -\nabla f_k$. Si ya diste pasos ($j>0$), sales usando el vector acumulado $p_k = z_j$ (que matemáticamente **garantiza** ser dirección de descenso).

---

## 6. Tabla «Chuleta» de 10 Segundos

| Si te preguntan por... | Responde esto |
| :--- | :--- |
| **Descenso de Gradiente** | Barato ($O(n)$), pero miope: zig-zag en valles estrechos, convergencia lineal lenta. |
| **Newton Exacto** | Velocísimo (cuadrático), pero carísimo ($O(n^3)$) y se rompe en puntos silla (Hessiana indefinida). |
| **DFP vs. BFGS** | Ambos son Cuasi-Newton de rango 2 y preservan $H_k \succ 0$, pero **BFGS se autocorrige** ante autovalores grandes y tolera búsquedas de línea inexactas (Wolfe), donde DFP suele frenarse. |
| **BFGS** | Solo usa gradientes, actualiza $H_k$ en $O(n^2)$ gracias a Sherman–Morrison–Woodbury, convergencia **superlineal**. Ideal para $n$ pequeño/mediano ($n \lesssim 10^3$). |
| **L-BFGS** | BFGS «a dieta»: en vez de guardar la matriz $H_k$ ($n \times n$), guarda solo los últimos $m \approx 10$ pares $(s_i, y_i)$. Memoria y costo $O(mn)$. Rey en Machine Learning clásico. |
| **Newton-CG** | Newton inexacto con Gradiente Conjugado. Memoria $O(n)$, no forma la Hessiana (*Hessian-Free*), evita *oversolving* con $\eta_k$ y trunca si $d_j^T H d_j \le 0$. Rey en gran escala con productos $H \cdot d$ baratos. |

---

## 7. Preguntas Trampa del Docente (y Cómo Responder Como un Crack)

### P1: *«¿Por qué BFGS actualiza $H_k$ (la inversa) y no $B_k$ directamente?»*
> **Respuesta:** Porque si actualizáramos $B_k$, en cada iteración tendríamos que resolver el sistema lineal $B_k p_k = -\nabla f_k$, lo cual cuesta $O(n^3)$ (igual que Newton). Al aplicar la fórmula de Sherman–Morrison–Woodbury actualizamos $H_k = B_k^{-1}$ directamente, y calcular $p_k = -H_k \nabla f_k$ baja a solo **$O(n^2)$**.

### P2: *«¿Qué pasa en BFGS si la condición de curvatura $s_k^T y_k > 0$ NO se cumple?»*
> **Respuesta:** Si $s_k^T y_k \le 0$, no existe ninguna matriz simétrica definida positiva que cumpla $B_{k+1} s_k = y_k$; la actualización BFGS perdería definitud positiva (o dividiría entre cero en $\rho_k = 1/(y_k^T s_k)$) y la dirección $p_k$ podría dejar de ser de descenso. Por eso usamos búsqueda de línea con **Condiciones de Wolfe**, que fuerzan $s_k^T y_k > 0$ automáticamente. (Como salvaguarda extra en código, si $s_k^T y_k \le 10^{-16}$ se omite la actualización en ese paso).

### P3: *«¿Por qué el método se llama BFGS y por qué desplazó a DFP?»*
> **Respuesta:** Por los cuatro autores que lo descubrieron de forma independiente en 1970: **Broyden, Fletcher, Goldfarb y Shanno**. Desplazó a DFP porque BFGS tiene una propiedad de **autocorrección**: si $H_k$ sobrestima mal la escala de un paso, BFGS corrige ese error en pocas iteraciones incluso con búsquedas de línea inexactas ($c_2 = 0.9$), mientras que DFP puede estancarse durante cientos de pasos (como se comprobó en nuestro experimento sobre Rosenbrock).

### P4: *«¿Cuál es la diferencia exacta entre Newton-Krylov (su primera expo) y Newton-CG?»*
> **Respuesta:** Newton-CG es un **caso particular** de Newton-Krylov especializado para **optimización irrestricta**. En sistemas generales $F(x)=0$ el Jacobiano no es simétrico y se usa GMRES. En optimización queremos $\nabla f(x)=0$, cuyo Jacobiano es la Hessiana $\nabla^2 f(x)$ que es **simétrica**, lo que permite usar **Gradiente Conjugado (CG)** (más barato en memoria y sin reinicios) añadiendo el **truncamiento por curvatura negativa** para garantizar que siempre descendamos hacia un mínimo y no hacia un punto silla.

### P5: *«¿Por qué truncar CG en $z_j$ cuando $d_j^T H d_j \le 0$ sigue siendo una dirección de descenso?»*
> **Respuesta:** Porque CG arranca en $z_0 = 0$ y en cada iteración anterior ($0, 1, \dots, j-1$) se verificó $d_i^T H d_i > 0$, lo que redujo estrictamente el modelo cuadrático $q(z) = \nabla f_k^T z + \frac{1}{2} z^T H_k z$. Como $q(0) = 0$, tenemos $q(z_j) < 0$; además $z_j^T H_k z_j > 0$ en el subespacio recorrido, lo que implica necesariamente $\nabla f_k^T z_j < 0$. Y si el fallo ocurre en $j=0$, tomamos $-\nabla f_k$, que trivialmente cumple $\nabla f_k^T (-\nabla f_k) = -\|\nabla f_k\|^2 < 0$.

### P6: *«Si tengo un problema con $n = 500{,}000$ variables, ¿uso BFGS o Newton-CG?»*
> **Respuesta:** Nunca BFGS denso (su matriz $H_k$ de $500{,}000 \times 500{,}000$ en `float64` necesitaría **2 Terabytes de RAM**). Las dos opciones reales son **Newton-CG** (solo $O(n) \approx 16$ MB de RAM, ideal si el producto $H \cdot d$ es barato o el problema está muy mal condicionado) o **L-BFGS** (memoria $O(mn) \approx 80$ MB con $m=10$, ideal si solo tengo gradientes).

### P7: *«¿Qué es el oversolving en Newton-CG y cómo se evita?»*
> **Respuesta:** Es gastar decenas o cientos de iteraciones internas de CG para resolver el sistema lineal con 12 dígitos de precisión cuando todavía estamos lejos de $x^*$ (donde la aproximación cuadrática de Taylor es mala). Se evita con una **secuencia de forzamiento adaptativa** como $\eta_k = \min(0.5, \sqrt{\|\nabla f_k\|})$, que pide poca precisión al principio y precisión máxima solo cuando ya estamos cerca del óptimo.

---

## 8. Guion Sugerido para una Exposición de 12–15 Minutos (2 Expositores)

| Tiempo | Expositor | Diapositivas | Qué enfatizar |
| :---: | :---: | :---: | :--- |
| **0:00 – 2:30** | Expositor 1 | Slides 1 – 3 | Conectar con Newton-Krylov; plantear el dilema GD vs. Newton Exacto; intuición geométrica de $H_k$ «redondeando» elipses. |
| **2:30 – 6:00** | Expositor 1 | Slides 4 – 7 | Ecuación secante $B_{k+1} s_k = y_k$, por qué Wolfe asegura $s_k^T y_k > 0$, evolución SR1 $\to$ DFP $\to$ BFGS, fórmula $H_{k+1}$ ($O(n^2)$ por Sherman–Morrison–Woodbury), convergencia superlineal y mención a L-BFGS. |
| **6:00 – 10:30** | Expositor 2 | Slides 8 – 11 | **Mover el slider interactivo** en Slide 8 ($n=10^6$); introducir Newton-CG, secuencia $\eta_k$ contra el *oversolving*, truco *Hessian-Free* y **ejecutar la demo en vivo del punto silla** en Slide 10; repasar el doble bucle de Newton-CG. |
| **10:30 – 13:30** | Ambos | Slides 12 – 14 | **Lanzar la carrera en vivo sobre Rosenbrock** en Slide 12 (mostrar los resultados del notebook); cerrar con la tabla comparativa, guía de decisión y conclusiones. |
