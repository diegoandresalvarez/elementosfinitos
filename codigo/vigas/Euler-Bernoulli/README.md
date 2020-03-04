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
Ver [aquí](../../repaso_matricial/portico_2d/deduccion_K_y_fe_elemento_portico_2D/).