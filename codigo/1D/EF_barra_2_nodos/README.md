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

## Deducción de las matriz de rigidez K y el vector de fuerzas nodales equivalentes
La salida de 
* [deduccion_K_y_f_EF_barra_2_nodos.m](deduccion_K_y_f_EF_barra_2_nodos.m)
* [deduccion_K_y_f_EF_barra_2_nodos.py](deduccion_K_y_f_EF_barra_2_nodos.py)

es:
```
K = 
  +-              -+
  |   E A     E A  |
  |   ---,  - ---  |
  |    L       L   |
  |                |
  |    E A   E A   |
  |  - ---,  ---   |
  |     L     L    |
  +-              -+

f = 
  +-     -+
  |  L b  |
  |  ---  |
  |   2   |
  |       |
  |  L b  |
  |  ---  |
  |   2   |
  +-     -+
```

## Ejemplo
Considere la barra mostrada a continuación:

![barra_con_carga_axial.svg](barra_con_carga_axial.svg)

### Solución mediante el método de los elementos finitos
* MATLAB: [barra_con_carga_axial_exacta_vs_EFs.m](barra_con_carga_axial_exacta_vs_EFs.m)
* PYTHON: [barra_con_carga_axial_exacta_vs_EFs.py](barra_con_carga_axial_exacta_vs_EFs.py)

### Solución resolviendo la ecuación diferencial asociada
NOTA: El siguiente programa hace uso de la función `bvp4c` de MATLAB.
* MATLAB: [barra_con_carga_axial_exacta_vs_bvp4c.m](barra_con_carga_axial_exacta_vs_bvp4c.m)
* PYTHON: [barra_con_carga_axial_exacta_vs_solve_bvp.py](barra_con_carga_axial_exacta_vs_solve_bvp.py)
