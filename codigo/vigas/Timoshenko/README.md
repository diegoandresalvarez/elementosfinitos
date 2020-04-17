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

El vector de fuerzas nodales equivalentes `fe` asociados a una carga distribuída `q` constante:
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
* [sec_2_4_EB_vs_T.m])(sec_2_4_EB_vs_T.m)

## Cálculo de las funciones de forma del elemento de viga de Euler-Bernoulli cuadrático
* [Kb_Ks_timoshenko_cuadratico.m](Kb_Ks_timoshenko_cuadratico.m)

## Ejemplo viga Euler-Bernoulli vs viga Timoshenko
eb_vs_t_h_20_1.png
c4_ejemplo_EB.m -> ../Euler-Bernoulli/c4_ejemplo_EB.m
c4_ejemplo_T.m


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
