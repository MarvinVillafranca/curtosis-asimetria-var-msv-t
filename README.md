# curtosis-asimetria-var-msv-t

Este repositorio contiene la implementación en:

### C++ (vía RcppArmadillo)
1. Funciones auxiliares `vec` e `inv_vec`.
2. Cálculo de matrices: \( V_t \), \( V_t^{-1} \), \( m_t \), \( M_t \), y \( S_t \).
3. Generación de la matriz de covarianza a priori de Minnesota.
4. Distribución condicional completa para estimar \( \boldsymbol{v}, A_1, \ldots, A_k \).
5. Función que calcula los residuales dados los parámetros estimados.
6. Cálculo de la matriz de covarianza de \( w_t \).
7. Cálculo de cuartos momentos (propiedades 4–8).
8. Cálculo de la curtosis de Koziol.
9. Matriz de covarianza de \( y_t \).

### R
1. Función que genera los parámetros del modelo.
2. Generación de los datos simulados.
   

## Archivos incluidos

- `momentos_GitHub.cpp`: Código en C++ que calcula la curtosis de Koziol y la asimetría de Mardia a partir de los parámetros estimados del modelo, usando `RcppArmadillo`.
- `simulacion_GitHub.R`: Código en R que contiene el proceso donde se generan los datos simulados.
 

## Descripción del proyecto

Este proyecto forma parte de una investigación sobre las propiedades de los cuartos momentos multivariados de los errores en 
modelos autorregresivos vectoriales con volatilidad estocástica multivariada y errores pesados (VAR-MSV-t). 
Se analiza cómo, bajo simulación, la asimetría de Mardia y la curtosis de Koziol convergen hacia sus valores teóricos.


## Requisitos

Este repositorio contiene funciones en C++ integradas con R mediante `RcppArmadillo`. Para compilar y ejecutar el código, necesitas:

- R (versión 4.0 o superior recomendada)
- Paquete `Rcpp`
- Paquete `RcppArmadillo`

Además, se asume que ya dispones de las matrices y vectores estimados de los parámetros latentes \(\alpha\), \(\Sigma\), \(\phi\), \(\nu\), los cuales deben estar disponibles antes de ejecutar las funciones contenidas en `momentos_GitHub.cpp`.

## Ejecución

Desde R, puedes compilar el archivo C++ de la siguiente manera:

```r
library(Rcpp)
sourceCpp("momentos_GitHub.cpp")

## Consideraciones sobre el código

Este repositorio implementa únicamente los aportes originales de los autores. Los parámetros latentes \(\alpha\), \(\Sigma\), \(\phi\) y \(\nu\), así como su proceso de estimación, provienen del artículo:

> Ishihara, T., & Omori, Y. (2012). Efficient Bayesian estimation of a multivariate stochastic volatility model with cross leverage and heavy-tailed errors. *Computational Statistics and Data Analysis*, 56(11), 3674–3693. https://doi.org/10.1016/j.csda.2010.07.015

Estos parámetros **no** se incluyen en este repositorio, ya que no forman parte de los desarrollos propios de esta investigación (aunque se propone dicha metodología como base para estimarlos).  
Se asume que están disponibles como entrada para las funciones aquí implementadas. Por tanto, **las funciones del archivo `momentos_GitHub.cpp` deben recibir estos valores como insumo** para el cálculo de los momentos de los errores y otras cantidades de interés.

## Publicación relacionada
Este repositorio acompaña el artículo:

"MODELO VAR INTEGRADO CON VOLATILIDAD ESTOCÁSTICA MULTIVARIADA Y ERRORES DE COLA PESADA",
Revista de Matemática Aplicada, Universidad de Costa Rica (en proceso de publicación).

El código aquí presentado corresponde a los desarrollos propios implementados como parte del artículo. Para reproducir completamente los resultados, se deben estimar los parámetros latentes siguiendo la metodología propuesta por Ishihara y Omori (2012), tal como se discute en el artículo.

