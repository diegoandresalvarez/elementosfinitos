# Programa para calcular los esfuerzos, deformaciones y desplazamientos de un sólido bidimensional utilizando elementos finitos triangulares de tres nodos

Calcule los campos de esfuerzos, desplazamientos y deformaciones de la viga mostrada:

![c5_viga_ejemplo.svg](c5_viga_ejemplo.svg)

Asuma:
* densidad del material = 7.8 kg/m³
* módulo de elasticidad = 200GPa
* coeficiente de Poisson = 0.30
* espesor de la viga = 10 cm

Se implementaron en PYTHON dos versiones, uno con matrices completas y otra con matrices ralas.

TIEMPO VERSION FULL:   0.0057590007781982 segundos (mejor tiempo de 5 corridas)
TIEMPO VERSION SPARSE: 0.0761265754699707 segundos (mejor tiempo de 5 corridas)
Es decir, la versión sparse() es 13 veces más lenta.

De otro lado
TAMAÑO MEMORIA "K" VERSION FULL:   749200 bytes
TAMAÑO MEMORIA "K" VERSION SPARSE:  41140 bytes
Es decir, para este ejemplo, la versión sparse() usa 18 veces menos memoria que la versión FULL.

En conclusión, solo se justifica utilizar SPARSE si estamos cortos de memoria                                                                                                                                                                                                             
