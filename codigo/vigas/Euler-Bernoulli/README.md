# Elementos finitos para modelar la flexión de vigas de Euler-Bernoulli

## Programas para el cálculo de las funciones de forma del elemento de viga de Euler-Bernoulli

### Usando polinomios interpoladores de Hermite
* MATLAB: [hermite.m](hermite.m)

### Resolviendo un sistema de ecuaciones
* MATLAB: [func_forma_euler_bernoulli.m](func_forma_euler_bernoulli.m). Nota: este programa adicionalmente calcula la matriz de rigidez `K` y la matriz de masa consistente `M`.

Estos programas verifican que las funciones de forma buscadas son:
```
         3
       xi    3 xi   1
N1 =   --- - ---- + -
        4      4    2

         3     2
       xi    xi    xi   1
N1b =  --- - --- - -- + -
        4     4     4   4

           3
         xi    3 xi   1
N2 =   - --- + ---- + -
          4      4    2

         3     2
       xi    xi    xi   1
N2b =  --- + --- - -- - -
        4     4     4   4
```
<img src="figs/funciones_forma_hermitianas.svg">

y que la matriz de rigidez del EF de viga de dos nodos es:
```
        /  12,  6 L,  -12,  6 L \
        |                       |
        |         2           2 |
    E I | 6 L, 4 L , -6 L, 2 L  |
K = ---*|                       |
      3 | -12, -6 L,  12,  -6 L |
     L  |                       |
        |         2           2 |
        \ 6 L, 2 L , -6 L, 4 L  /
```

## Un polinomio de orden `n` y otro de orden `n-1` ajustado por mínimos cuadrados se intersectan en la raíces del polinomio de Legendre de orden `n`

Cuando la intersección de un polinomio de grado `n` con su ajuste por mínimos cuadrados un polinomio de grado `n-1` sucede, se pueden observar que dicha intersección ocurre en raíces del polinomio de Legendre de orden `n`:

<img src="../../2D/extrapolacion_de_esfuerzos/figs/interseccion_polinomios_en_raices_pol_Legendre.png">

Código:
* MATLAB: [interseccion_polinomios_en_raices_pol_Legendre.m](../../2D/extrapolacion_de_esfuerzos/interseccion_polinomios_en_raices_pol_Legendre.m)

##  Cálculo de la matriz K para el EF de 2 nodos de Euler-Bernoulli resolviendo la ecuación diferencial
* MATLAB: [K_EB_exacta.m](K_EB_exacta.m)
Ver [aquí](../../repaso_matricial/portico_2d/deduccion_K_y_fe_elemento_portico_2D/).

## Vector f asociado a una carga trapezoidal
Con el siguiente código se estableció que las fuerzas nodales equivalentes son:
* [f_EB_carga_trapezoidal.m](f_EB_carga_trapezoidal.m)
<img src="figs/viga_carga_trapezoidal.svg">

## Vector f asociado a los casos 3 y 6 (diapositivas de Logan)
Como aquí se requieren funciones de discontinuidad para resolver la ecuación diferencial, se utilizó MAXIMA, ya que MATLAB tiene un bug con funciones de discontinuidad y en SYMPY estos casos aun no están implementados. Recuerde que el código usa la librería [pw](https://sourceforge.net/projects/piecewisefunc/).
* [f_casos_3_y_6.mac](f_casos_3_y_6.mac)

## Cálculo de los diagramas de cortante, momento, ángulo de inclinación y desplazamiento de una viga

Considere la siguiente viga, la cual se tomó del libro [J. Uribe Escamilla (2002) - Análisis de estructuras](https://www.researchgate.net/publication/31754481_Analisis_de_estructuras_J_Uribe_Escamilla):

<img src="resolviendo_la_ecuacion_diferencial/c4_escamilla_ej_5_5.png">

Dicha viga se resolvió utilizando la teoría de **Euler-Bernoulli**:

* Solución resolviendo directamente la ecuación diferencial (usando la funcion bvp5c): 
  * MATLAB: [c4_escamilla_ej_5_5_EB_eq_diff.m]

* Solución por el método de los elementos finitos:
  * MATLAB: [c4_ejemplo_EB.m](c4_ejemplo_EB.m)

La salida de este último código es:

<img src="../ejemplos/figs/c4_ejemplo_EB_v_theta.svg">

<img src="../ejemplos/figs/c4_ejemplo_EB_M_V.svg">

siendo las reacciones de la viga:
```
Fuerzas nodales de equilibrio (solo se imprimen los diferentes de cero)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Nodo   1 Ry =         15.4 kN, Mz =        38.83 kN-m
Nodo 101 Ry =        48.91 kN, Mz =            0 kN-m
Nodo 161 Ry =        52.69 kN, Mz =            0 kN-m
```
## Cálculo de vigas con rótulas

[viga_con_rotulas/](viga_con_rotulas)