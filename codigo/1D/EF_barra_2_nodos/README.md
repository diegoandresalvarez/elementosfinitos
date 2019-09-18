# Elemento finito de barra de dos nodos

El elemento finito considerado es el siguiente:

![EF_barra_2_nodos.svg](EF_barra_2_nodos.svg)

se asume que está sometido a una carga axial de magnitud constante.

## Deducción de las funciones de forma

La salida de 
* [deduccion_func_forma_EF_barra_2_nodos.m](deduccion_func_forma_EF_barra_2_nodos.m)
* [deduccion_func_forma_EF_barra_2_nodos.py](deduccion_func_forma_EF_barra_2_nodos.py)

es:
```
          u1 x2 - u2 x1
a0 =   - -------------
            x1 - x2

       u1 - u2
a1 =   -------
       x1 - x2

      /    x        x2    \      /   x1         x    \
u =   | ------- - ------- | u1 + | ------- - ------- | u2
      \ x1 - x2   x1 - x2 /      \ x1 - x2   x1 - x2 /
```

## Deducción de la matriz de rigidez `K` y el vector de fuerzas nodales equivalentes `f`
La salida de 
* [deduccion_K_y_f_EF_barra_2_nodos.m](deduccion_K_y_f_EF_barra_2_nodos.m)
* [deduccion_K_y_f_EF_barra_2_nodos.py](deduccion_K_y_f_EF_barra_2_nodos.py)

es:
```
K = 
/  A E     A E \
|  ---,  - --- |
|   L       L  |
|              |
|   A E   A E  |
| - ---,  ---  |
\    L     L   /

f = 
/ L b \
| --- |
|  2  |
|     |
| L b |
| --- |
\  2  /
```

## Ejemplo
Considere la barra mostrada a continuación:

![barra_con_carga_axial.svg](barra_con_carga_axial.svg)

### Solución teórica utilizando funciones de forma globales
* MATLAB: [deduccion_K_y_f_func_forma_global.m](deduccion_K_y_f_func_forma_global.m)
* PYTHON: [deduccion_K_y_f_func_forma_global.py](deduccion_K_y_f_func_forma_global.py)

La salida del programa es la siguiente:
```
K = 
    ⎡   L                                                                                                  ⎤
    ⎢   ⌠                      L                         L                         L                       ⎥
    ⎢   ⎮            2         ⌠                         ⌠                         ⌠                       ⎥
    ⎢   ⎮ ⎛d        ⎞          ⎮ d         d             ⎮ d         d             ⎮ d         d           ⎥
    ⎢   ⎮ ⎜──(N₁(x))⎟  dx      ⎮ ──(N₁(x))⋅──(N₂(x)) dx  ⎮ ──(N₁(x))⋅──(N₃(x)) dx  ⎮ ──(N₁(x))⋅──(N₄(x)) dx⎥
    ⎢   ⎮ ⎝dx       ⎠          ⎮ dx        dx            ⎮ dx        dx            ⎮ dx        dx          ⎥
    ⎢   ⌡                      ⌡                         ⌡                         ⌡                       ⎥
    ⎢   0                      0                         0                         0                       ⎥
    ⎢                                                                                                      ⎥
    ⎢                             L                                                                        ⎥
    ⎢L                            ⌠                      L                         L                       ⎥
    ⎢⌠                            ⎮            2         ⌠                         ⌠                       ⎥
    ⎢⎮ d         d                ⎮ ⎛d        ⎞          ⎮ d         d             ⎮ d         d           ⎥
    ⎢⎮ ──(N₁(x))⋅──(N₂(x)) dx     ⎮ ⎜──(N₂(x))⎟  dx      ⎮ ──(N₂(x))⋅──(N₃(x)) dx  ⎮ ──(N₂(x))⋅──(N₄(x)) dx⎥
    ⎢⎮ dx        dx               ⎮ ⎝dx       ⎠          ⎮ dx        dx            ⎮ dx        dx          ⎥
    ⎢⌡                            ⌡                      ⌡                         ⌡                       ⎥
    ⎢0                            0                      0                         0                       ⎥
A⋅E⋅⎢                                                                                                      ⎥
    ⎢                                                       L                                              ⎥
    ⎢L                         L                            ⌠                      L                       ⎥
    ⎢⌠                         ⌠                            ⎮            2         ⌠                       ⎥
    ⎢⎮ d         d             ⎮ d         d                ⎮ ⎛d        ⎞          ⎮ d         d           ⎥
    ⎢⎮ ──(N₁(x))⋅──(N₃(x)) dx  ⎮ ──(N₂(x))⋅──(N₃(x)) dx     ⎮ ⎜──(N₃(x))⎟  dx      ⎮ ──(N₃(x))⋅──(N₄(x)) dx⎥
    ⎢⎮ dx        dx            ⎮ dx        dx               ⎮ ⎝dx       ⎠          ⎮ dx        dx          ⎥
    ⎢⌡                         ⌡                            ⌡                      ⌡                       ⎥
    ⎢0                         0                            0                      0                       ⎥
    ⎢                                                                                                      ⎥
    ⎢                                                                                 L                    ⎥
    ⎢L                         L                         L                            ⌠                    ⎥
    ⎢⌠                         ⌠                         ⌠                            ⎮            2       ⎥
    ⎢⎮ d         d             ⎮ d         d             ⎮ d         d                ⎮ ⎛d        ⎞        ⎥
    ⎢⎮ ──(N₁(x))⋅──(N₄(x)) dx  ⎮ ──(N₂(x))⋅──(N₄(x)) dx  ⎮ ──(N₃(x))⋅──(N₄(x)) dx     ⎮ ⎜──(N₄(x))⎟  dx    ⎥
    ⎢⎮ dx        dx            ⎮ dx        dx            ⎮ dx        dx               ⎮ ⎝dx       ⎠        ⎥
    ⎢⌡                         ⌡                         ⌡                            ⌡                    ⎥
    ⎣0                         0                         0                            0                    ⎦

f = 
  ⎡L         ⎤
  ⎢⌠         ⎥
  ⎢⎮ N₁(x) dx⎥
  ⎢⌡         ⎥
  ⎢0         ⎥
  ⎢          ⎥
  ⎢L         ⎥
  ⎢⌠         ⎥
  ⎢⎮ N₂(x) dx⎥
  ⎢⌡         ⎥
  ⎢0         ⎥
b⋅⎢          ⎥
  ⎢L         ⎥
  ⎢⌠         ⎥
  ⎢⎮ N₃(x) dx⎥
  ⎢⌡         ⎥
  ⎢0         ⎥
  ⎢          ⎥
  ⎢L         ⎥
  ⎢⌠         ⎥
  ⎢⎮ N₄(x) dx⎥
  ⎢⌡         ⎥
  ⎣0         ⎦

```

### Solución mediante el método de los elementos finitos
* MATLAB: [barra_con_carga_axial_exacta_vs_EFs.m](barra_con_carga_axial_exacta_vs_EFs.m)
* PYTHON: [barra_con_carga_axial_exacta_vs_EFs.py](barra_con_carga_axial_exacta_vs_EFs.py)

### Solución resolviendo la ecuación diferencial asociada
* MATLAB: [barra_con_carga_axial_exacta_vs_bvp4c.m](barra_con_carga_axial_exacta_vs_bvp4c.m)
* PYTHON: [barra_con_carga_axial_exacta_vs_solve_bvp.py](barra_con_carga_axial_exacta_vs_solve_bvp.py)
