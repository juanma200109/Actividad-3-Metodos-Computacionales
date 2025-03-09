# Actividad-3-Metodos-Computacionales
# Flujo de Potencia Cuasi Dinámico usando el Método de Punto Fijo
*Este proyecto implementa un análisis de flujo de potencia cuasi dinámico en sistemas eléctricos utilizando el método de punto fijo. El objetivo es determinar el estado del sistema, específicamente las tensiones en los nodos, considerando variaciones temporales en las demandas y la generación. El algoritmo resuelve las ecuaciones no lineales del flujo de potencia de manera iterativa y genera resultados que pueden usarse para análisis posteriores, como la estabilidad del sistema o la planificación operativa.*

----

## Marco Teórico

*El flujo de potencia es un análisis esencial en sistemas eléctricos que permite encontrar el estado del sistema (tensiones y flujos de potencia) a partir de valores nodales. Este proyecto utiliza el método de punto fijo para resolver las ecuaciones no lineales asociadas.*

#### Inputs
- **Datos de líneas**: Resistencias, reactancias, susceptancias y taps de transformadores.
- **Datos de nodos**: Tipo de nodo (Slack, PV, PQ), potencias generadas ($P_g$, $Q_g$) y consumidas ($P_d$, $Q_d$).
- **Demandas**: Perfiles de carga diarios o horarios para cada nodo tipo PQ.
- **Generación**: Perfiles de generación diarios o horarios, como energía solar.

#### Outputs
- **Tensiones en nodos**: Calculadas para cada nodo no-Slack en cada intervalo de tiempo.
- Archivos CSV con las tensiones para análisis posteriores.

#### Clasificación de Nodos
- **Nodos Slack**: Representa la subestación o barra fija (tensión controlada, usualmente 1∠0°). Absorbe las diferencias de potencia del sistema.
- **Nodos PV**: Generadores que controlan su tensión y potencia activa.
- **Nodos PQ**: Nodos donde su potencia activa y reactiva son conocidas y controladas.

#### Método de Punto Fijo
*El método de punto fijo resuelve iterativamente las ecuaciones del flujo de potencia en un espacio de Banach completo (e.g., $\mathbb{C}^n, \mathbb{R}^n$). Se define una transformación $T: \mathbb{C}^n \rightarrow \mathbb{C}^n$ como un mapa de contracción si satisface:
$||T(x)-T(y)||<k||x-y||$, donde $k<1$ es una constante de contracción.

*Esto garantiza un único punto fijo, calculado mediante:*
$x^{(k+1)} = T(x^{(k)})$

En el contexto del flujo de potencia:
1. Los nodos se dividen en:
    - **Nodos S**: Nodo Slack (subestación).
    - **Nodos N**: Nodos restantes (PV y PQ).
2. La ecuación base es: 
    $$
\begin{bmatrix}
\textit{I}_s \\
\textit{I}_n
\end{bmatrix}
=
\begin{bmatrix}
\textit{Y}_{ss} & \textit{Y}_{sn} \\
\textit{Y}_{ns} & \textit{Y}_{nn}
\end{bmatrix}
\begin{bmatrix}
\textit{V}_s \\
\textit{V}_n
\end{bmatrix}
$$

3. Las potencias nodales se relacionan con las tensiones:
$S_n=V_n \cdot I_n^*$

4. La ecuación iterativa para $V_n$ es: $$V_n^{k+1}=Y_{nn}^{-1}\left(\left(\frac{S_n}{V_n^{(k)}}\right)^*-Y_{ns}V_s\right)$$

    *Donde*:
    - $Y_{nn}$ es la matriz de admisión de los nodos no-Slack.
    - $Y_{ns}$ es la matriz de admisión de los nodos no-Slack a los nodos Slack.
    - $S_n$: Vector de potencias netas inyectadas.
    - $V_n$: Vector de tensiones de los nodos no-Slack.
    - $V_s$: Tensión del nodo Slack.

#### Criterio de Convergencia
El proceso iterativo se detiene cuando la norma del error entre dos iteraciones consecutivas es menor que un umbral (e.g., $1*10^{-6}$):
$||V_n^{k+1}-V_n^{k}||<1*10^{-6}$
La norma se define como:
$||V_n|| = \sqrt{V_1^2 + V_2^2 + \cdots + V_n^2}$

---
## **Funciones:**
*Para ejecutar este código, necesitas lo siguiente:*
- *Julia 1.11 o superior*
- *Bibliotecas*:
    - LinearAlgebra
    - Clusterinig
    - Statistics
    - DataFrames
    - Plots
    - CSV

**calcular_ynn_yns:**

*Descripción*

Construye la matriz de admitancias $(Y_{bus})$ del sistema de potencia y extrae las submatrices relevantes para estudios de flujo de potencia. Considera parámetros de líneas como resistencia, reactancia, susceptancia shunt y ajustes de tap en transformadores. La matriz Ybus es esencial para análisis de flujo de potencia y estabilidad.

*Requiere*

    - using LinearAlgebra
    - using DataFrames

    Entradas:  

    - lines: DataFrame con la información de las líneas del sistema. Contiene columnas como:
    
      - FROM: Nodo de envío (ej. "N1").
      - TO: Nodo de recibo (ej. "N2").
      - R: Resistencia de la línea.
      - x: Resistencia de la línea.
      - b: Susceptancia shunt de la línea.
      - tap: Ajuste de tap en transformadores.

    - nodes : DataFrame con la información de los nodos. Contiene columnas como:

      - NUMBER: Número del nodo (ej. "N1").
      - TYPE: Tipo de nodo (3 indica el nodo slack).

    Salida :    

    - Y_nn: Matriz de admitancias (compleja) de los nodos excluyendo el slack.
    - Y_ns: Vector de admitancias (complejo) asociado al nodo slack.
    - Ybus: Matriz de admitancias completa del sistema.

**flujo_punto_fijo:**

*Descripción*

Calcula las tensiones en los nodos mediante el método iterativo de punto fijo. Utiliza la matriz $Y_{nn}$ y el vector $Y_{ns}$ para resolver el sistema de ecuaciones no lineales. Incluye un gráfico de convergencia del error en escala logarítmica.

*Requiere*

    - using LinearAlgebra
    - using DataFrames
    - using Plots

    Entradas:  

    - lines: DataFrame con la información de las líneas del sistema. Contiene columnas como:
    
      - FROM: Nodo de envío (ej. "N1").
      - TO: Nodo de recibo (ej. "N2").
      - R: Resistencia de la línea.
      - x: Resistencia de la línea.
      - b: Susceptancia shunt de la línea.
      - tap: Ajuste de tap en transformadores.

    - nodes : DataFrame con la información de los nodos. Contiene columnas como:

      - NUMBER: Número del nodo (ej. "N1").
      - TYPE: Tipo de nodo (3 indica el nodo slack).

    - Bbus  : Matriz de susceptancia reducida, obtenida de B_bus.

    Salida :    

    - Vn: Vector de tensiones complejas en los nodos no slack.
    - Gráfico de convergencia del error (iteraciones vs. error en escala logarítmica).

**flujo_cuasi_dinamico:**

*Descripción*

Simula un flujo de potencia cuasi-dinámico para múltiples escenarios (días típicos) e intervalos de tiempo de 5 minutos. Calcula las tensiones nodales considerando variaciones en generación y demanda, y almacena los resultados en archivos CSV.

*Requiere*

    - using LinearAlgebra
    - using DataFrames
    - using CSV

    Entradas:  

    - lines: DataFrame con la información de las líneas del sistema. Contiene columnas como:
    
      - FROM: Nodo de envío (ej. "N1").
      - TO: Nodo de recibo (ej. "N2").
      - R: Resistencia de la línea.
      - x: Resistencia de la línea.
      - b: Susceptancia shunt de la línea.
      - tap: Ajuste de tap en transformadores.

    - nodes : DataFrame con la información de los nodos. Contiene columnas como:

      - NUMBER: Número del nodo (ej. "N1").
      - TYPE: Tipo de nodo (3 indica el nodo slack).

    - Demandas_prom: DataFrame con demandas promedio por nodo y tiempo.
    - dias_tipicos: DataFrame con perfiles de generación para días típicos (columnas = casos, filas = intervalos).

    Salida :    

    - Archivos CSV (ej. Vn_caso_1.csv) con las tensiones calculadas para cada caso e intervalo de tiempo.


#### Notas:
- *Puedes instalar las bibliotecas en Julia ejecutando:*
    ```julia
    using Pkg
    Pkg.add(["LinearAlgebra", "DataFrames", "CSV", "Plots", "Dates", "Clustering", "Statistics"])
- Las funciones asumen que el nodo slack está identificado con TYPE = 3 en el DataFrame nodes.
- El método de punto fijo en flujo_punto_fijo se limita a 4 iteraciones para demostración.
- flujo_cuasi_dinamico utiliza un enfoque simplificado con iteraciones fijas (sin criterio de convergencia dinámico).

[![License: CC BY-NC-ND 4.0](https://img.shields.io/badge/License-CC_BY--NC--ND_4.0-lightgrey)](https://creativecommons.org/licenses/by-nc-nd/4.0/)
