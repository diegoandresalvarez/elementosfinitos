# Programa para calcular los esfuerzos, deformaciones y desplazamientos de un sólido bidimensional utilizando elementos finitos isoparamétricos rectangulares de ocho nodos

Calcule los campos de esfuerzos, desplazamientos y deformaciones de la estructura mostrada:
![c5_isoparametric_cuad_8_nodos.jpg](c5_isoparametric_cuad_8_nodos.jpg)

* MATLAB: [c5_ejemplo_isoparametricos_rect_8_nodos.zip] 
* JULIA 0.5.1 (experimental): [c5_ejemplo_isoparametricos_rect_8_nodos_julia_0.51.zip]

Los esfuerzos de sigma_1 y sigma_2 calculados por el programa son:

![c5_ejemplo_isoparametricos_rect_8_nodos_s1_s2.png](c5_ejemplo_isoparametricos_rect_8_nodos_s1_s2.png)

NOTA 1: la matriz para extrapolar los esfuerzos desde los puntos de Gauss hacia los nodos se dedujo con el programa 
* MATLAB: [c5_extrapolacion_esfuerzos.m]

NOTA 2: falta calcular las deformaciones principales. 

NOTA 3: este programa exporta los resultados a [GiD](http://gid.cimne.upc.es/] por ejemplo: 

![c5_ejemplo_isoparametricos_rect_8_nodos_exportar_resultados_gid.png](c5_ejemplo_isoparametricos_rect_8_nodos_exportar_resultados_gid.png)

Con el mismo programa, escogiendo la malla 3 (creada por Alejandro Cardona Jimenez), calculamos:

![c5_gancho.png](c5_gancho.png)

