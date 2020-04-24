# Elementos finitos para modelar la flexión de vigas de Timoshenko

## Cálculo de las funciones de forma del elemento de viga de Euler-Bernoulli lineal
El programa:
* [Kb_Ks_timoshenko_lineal.m](Kb_Ks_timoshenko_lineal.m)

calcula:

Las matrices de deformación por flexión `Bb` y por cortante `Bs`:

```
Bb =
[ 0, -1/L, 0, 1/L]
```
```
Bs =
[ -1/L, xi/2 - 1/2, 1/L, - xi/2 - 1/2]
```

Las matrices `Kb` y `Ks` que se calculan por la integración exacta:
```
Kb = (E*I/L) * 
/ 0,  0, 0,  0 \
|              |
| 0,  1, 0, -1 |
|              |
| 0,  0, 0,  0 |
|              |
\ 0, -1, 0,  1 /
```
```
Ks = (G*Aast/L) * 
/      L         L  \
|  1,  -,   -1,  -  |
|      2         2  |
|                   |
|       2         2 |
|  L   L     L   L  |
|  -,  --, - -,  -- |
|  2    3    2    6 |
|                   |
|       L         L |
| -1, - -,  1,  - - |
|       2         2 |
|                   |
|       2         2 |
|  L   L     L   L  |
|  -,  --, - -,  -- |
\  2    6    2    3 /
```

Las matrices `Kb` y `Ks` que se calculan con Gauss-Legendre usando 1 punto de integración:
```
Kb = (E*I/L) * 
/ 0,  0, 0,  0 \
|              |
| 0,  1, 0, -1 |
|              |
| 0,  0, 0,  0 |
|              |
\ 0, -1, 0,  1 /
```
```
Ks = (G*Aast/L) * 
/      L         L  \
|  1,  -,   -1,  -  |
|      2         2  |
|                   |
|       2         2 |
|  L   L     L   L  |
|  -,  --, - -,  -- |
|  2    4    2    4 |
|                   |
|       L         L |
| -1, - -,  1,  - - |
|       2         2 |
|                   |
|       2         2 |
|  L   L     L   L  |
|  -,  --, - -,  -- |
\  2    4    2    4 /
```

Las matrices `Kb` y `Ks` que se calculan con Gauss-Legendre usando 2 puntos de integración:
```
Kb = (E*I/L) * 
/ 0,  0, 0,  0 \
|              |
| 0,  1, 0, -1 |
|              |
| 0,  0, 0,  0 |
|              |
\ 0, -1, 0,  1 /
```
```
Ks = (G*Aast/L) * 
/      L         L  \
|  1,  -,   -1,  -  |
|      2         2  |
|                   |
|       2         2 |
|  L   L     L   L  |
|  -,  --, - -,  -- |
|  2    3    2    6 |
|                   |
|       L         L |
| -1, - -,  1,  - - |
|       2         2 |
|                   |
|       2         2 |
|  L   L     L   L  |
|  -,  --, - -,  -- |
\  2    6    2    3 /
```

El vector de fuerzas nodales equivalentes `fe` asociados a una carga distribuída `q` constante y un momento distribuído `m` tambien uniforme:
```
fe = 
/ L fz \
| ---- |
|   2  |
|      |
|  L m |
|  --- |
|   2  |
|      |
| L fz |
| ---- |
|   2  |
|      |
|  L m |
|  --- |
\   2  /
```
El vector de fuerzas nodales equivalentes `fe` asociados a una carga distribuída trapezoidal y sin momentos distribuidos:
```
fe = 
/ L (2 q1 + q2) \
| ------------- |
|       6       |
|               |
|       0       |
|               |
| L (q1 + 2 q2) |
| ------------- |
|       6       |
|               |
\       0       /
```

## Ejercicio de la Sección 2.4 de Oñate: comparación EB vs T (1 GL) vs T (2 GL)
El código

* [sec_2_4_EB_vs_T.m](sec_2_4_EB_vs_T.m)

implementa el ejercicio que hay en la Sección 2.4 de Oñate (2013) y obtiene este gráfico:

<img src="figs/lambda_vs_rw_0_a_6.svg">

Básicamente este ejercicio compara los desplazamientos calculados con las teorías de Euler-Bernoulli y Timoshenko utilizando una matriz de rigidez por cortante `Ks` calculada con integración reducida (1 punto de GL) e integración completa (2 puntos de GL). Se deduce de este ejercicio que:
* En una viga esbelta (`λ` muy grande) el efecto del esfuerzo cortante es despreciable y la solución numérica coincide con la predicha por la teoría de Euler-Bernoulli.
* A medida que el número de EFs aumenta, también así lo hace la calidad de la solución para el EF con que aplica integración reducida a `Ks`.
* Con la integración exacta de `Ks` se produce el fenómeno de *shear locking* (bloqueo de la solución). Dicho fenómeno hace que la viga sea en el límite `λ→∞` infinitamente rígida. Este EF con integración exacta solo funcionaría con un número exagerado de EFs (asumiendo que `λ` no tiende a infinito), y aún así su precisión no sería buena, lo que lo hace inutilizable en la mayoría de los casos.
* Con la integración reducida de `Ks` se evita el fenómeno del bloqueo por cortante (shear locking) y el EF resultante es válido para vigas de pequeño y gran canto. Veremos más adelante, que en este caso el punto central de Gauss-Legendre es adicionalmente el punto óptimo para el cálculo de los esfuerzos.

## Ejemplo viga Euler-Bernoulli vs viga Timoshenko
El programa 
* [c4_ejemplo_T.m](c4_ejemplo_T.m)
* [c4_ejemplo_EB.m -> ../Euler-Bernoulli/c4_ejemplo_EB.m](c4_ejemplo_EB.m)
calcula la viga de Timoshenko con elementos finitos de dos nodos. Dicho programa hace una comparación con el método de Euler-Bernoulli

Por ejemplo, el programa anterior calcula la viga:

<img src="../ejemplos/figs/viga_Uribe_Escamilla_ej_5_5.png">

obteniendo la siguiente comparación (para h = 2.0 m):

<img src="figs/eb_vs_t_h_20_1.png">


## Cálculo de las funciones de forma del elemento de viga de Euler-Bernoulli cuadrático
* [Kb_Ks_timoshenko_cuadratico.m](Kb_Ks_timoshenko_cuadratico.m)

##  Cálculo de la matriz K para el EF de 2 nodos calculado utilizando integración exacta:
* [K_exacta_viga_T.m](K_exacta_viga_T.m)
* [f_exacta_carga_trapezoidal_T.m](f_exacta_carga_trapezoidal_T.m)
* [c4_ejemplo_con_K_T_exacta.m](c4_ejemplo_con_K_T_exacta.m)

## Interpolación acoplada (linked interpolation):
* [ej_2_2_interpolacion_acoplada.m](ej_2_2_interpolacion_acoplada.m)
* [ej_2_3.m](ej_2_3.m)
* [ej_2_4.m](ej_2_4.m)
* [ej_2_5.m](ej_2_5.m)
* [ej_2_6.m](ej_2_6.m)
* [ej_2_7.m](ej_2_7.m)


## Imposición de un campo de deformaciones angulares para gxz:
* [Bs_sustitutiva_T2.m](Bs_sustitutiva_T2.m)
* [Bs_sustitutiva_T3.m](Bs_sustitutiva_T3.m)
* [ej_2_8_imposicion_gxz_T2.m](ej_2_8_imposicion_gxz_T2.m)
* [ej_2_9_imposicion_gxz_T3.m](ej_2_9_imposicion_gxz_T3.m)
