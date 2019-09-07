# Elemento finito de barra de dos nodos

El elemento finito considerado es el siguiente:

![c2_barra_con_carga.svg](c2_barra_con_carga.svg)

se asume que está sometido a una carga axial de magnitud constante.

## Deducción de las funciones de forma

La salida de [c2_deduccion_func_forma_EF_barra_2_nodos.m](c2_deduccion_func_forma_EF_barra_2_nodos.m) es:
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
La salida de [c2_deduccion_K_y_f_EF_barra_2_nodos.m](c2_deduccion_K_y_f_EF_barra_2_nodos.m) es:
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

![barraempotrada.svg](barraempotrada.svg)

### Solución mediante el método de los elementos finitos
[[image:c3_ejemplo_barra.png]]
* [c2_ejemplo_barra_con_carga_axial_exacta_vs_EFs.m](c2_ejemplo_barra_con_carga_axial_exacta_vs_EFs.m)

### Solución resolviendo la ecuación diferencial asociada
El siguiente programa hace uso de la función `bvp4c` de MATLAB.
* [c2_ejemplo_barra_con_carga_axial_exacta_vs_bvp4c.m](c2_ejemplo_barra_con_carga_axial_exacta_vs_bvp4c.m)
