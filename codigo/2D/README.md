**<span style="font-size: 300%">CODIGOS DE MATLAB</span>**
<span style="color: #0000ff;
font-size: 200%;">Nota: estos códigos están hechos para que funcionen con MATLAB 2013a. No los he actualizado a versiones más nuevas de MATLAB, ya que el 2013a es el MATLAB que se tiene instalado en los computadores de la Universidad Nacional de Colombia - Sede Manizales. Por lo tanto, los programas para álgebra simbólica podrían fallar si usted utiliza versiones modernas de MATLAB.</span>

[[image:http://imgs.xkcd.com/comics/ballmer_peak.png]]
Fuente: [[http://xkcd.com/323/]]

=CAPITULO 5=

==Deducción de las funciones de forma de un elemento triangular de tres nodos==
Con la ayuda de este programa de MATLAB: [[file:c5_deduccion_func_forma_triangulo_3_nodos.m]] se pueden encontrar dichas funciones.
La ejecución de dicho código se puede ver al dar click [[http://elementosfinitosunalmzl.wikispaces.com/file/view/c5_deduccion_func_forma_triangulo_3_nodos.html|AQUI]]. El resultado de la ejecución es:
[[code]]
u = 
   / x2 y3 - x3 y2   y (x2 - x3)   x (y2 - y3) \
   | ------------- - ----------- + ----------- | u1 +
   \    2 Area         2 Area        2 Area    /

   / y (x1 - x3)   x1 y3 - x3 y1   x (y1 - y3) \
   | ----------- - ------------- - ----------- | u2 +
   \   2 Area         2 Area         2 Area    /

   / x1 y2 - x2 y1   y (x1 - x2)   x (y1 - y2) \
   | ------------- - ----------- + ----------- | u3
   \    2 Area         2 Area        2 Area    /
[[code]]


==Programa para calcular los esfuerzos, deformaciones y desplazamientos de un sólido bidimensional utilizando elementos finitos triangulares de tres nodos==
Considere la viga mostrada, suponiendo que el peso del material es 7.8 kg/m3, E = 200GPa, el coeficiente de Poisson es 0.30 y el espesor de la viga es 10 cm. Calcule los campos de esfuerzos, desplazamientos y deformaciones de la viga
[[image:c5_viga_ejemplo.png]]
Código de MATLAB: [[file:c5_ejemplo_triang_3_nodos_k_dado.zip]] <span style="color: #ff0000;">(NOTA: el código está bastante particularizado al ejemplo. Debe tenerse cuidado si se quiere utilizar para otro tipo de estructura)</span>


==Deducción de la matriz de rigidez de un elemento rectangular de 4 nodos==
Con la ayuda de este programa de MATLAB: [[file:c5_deduccion_func_forma_rectangulo_4_nodos.m]] se pueden encontrar K.

El resultado de la ejecución es:
[[code]]
K =
[ 2*a1 + 2*a4,         a36,         c41,         b36,        -a14,        -a36,         c14,         b63]
[         a36, 2*a2 + 2*a5,         b63,         c25,        -a36,        -a25,         b36,         c52]
[         c41,         b63, 2*a1 + 2*a4,        -a36,         c14,         b36,        -a14,         a36]
[         b36,         c25,        -a36, 2*a2 + 2*a5,         b63,         c52,         a36,        -a25]
[        -a14,        -a36,         c14,         b63, 2*a1 + 2*a4,         a36,         c41,         b36]
[        -a36,        -a25,         b36,         c52,         a36, 2*a2 + 2*a5,         b63,         c25]
[         c14,         b36,        -a14,         a36,         c41,         b63, 2*a1 + 2*a4,        -a36]
[         b63,         c52,         a36,        -a25,         b36,         c25,        -a36, 2*a2 + 2*a5]

 
Para facilitar la comparacion de la matriz con el libro de Oñate
se muestra únicamente la parte superior de K =
[ 2*a1 + 2*a4,         a36,         c41,         b36,        -a14,        -a36,         c14,         b63]
[           0, 2*a2 + 2*a5,         b63,         c25,        -a36,        -a25,         b36,         c52]
[           0,           0, 2*a1 + 2*a4,        -a36,         c14,         b36,        -a14,         a36]
[           0,           0,           0, 2*a2 + 2*a5,         b63,         c52,         a36,        -a25]
[           0,           0,           0,           0, 2*a1 + 2*a4,         a36,         c41,         b36]
[           0,           0,           0,           0,           0, 2*a2 + 2*a5,         b63,         c25]
[           0,           0,           0,           0,           0,           0, 2*a1 + 2*a4,        -a36]
[           0,           0,           0,           0,           0,           0,           0, 2*a2 + 2*a5]
[[code]]


==Cálculo de las funciones de forma del elemento rectangular lagrangiano de 16 nodos ==
Código de MATLAB: [[file:c5_funciones_forma_lagrangianos_rect_2D.m]]

Con este código se obtuvo por ejemplo que la función de forma 12 de este elemento es:
[[code]]
N12 = 
    9 (xi - 1) (3 xi - 1) (3 xi + 1) (eta - 1) (3 eta - 1) (eta + 1)
  - ----------------------------------------------------------------
                                  256
[[code]]

Siendo el gráfico de esta función
[[image:c5_func_forma_N12_rect_16_nodos.png]]


==Cálculo de las funciones de forma del elemento rectangular serendípito de 8 nodos ==
Código de MATLAB: [[file:c5_funciones_forma_serendipito_rect_2D_8_nodos.m]]

Con este código se obtuvo por ejemplo que la función de forma 4 de este elemento es:
[[code]]
N4 = 
        2
    (eta  - 1) (xi + 1)
  - -------------------
             2
[[code]]

Siendo el gráfico de esta función:
[[image:c5_func_forma_N4_rect_serend_8_nodos.png]]


==Cálculo de las funciones de forma del elemento triangular de 10 nodos==
Código de MATLAB: [[file:c5_funciones_forma_triang.m]]

Con este código se obtuvo por ejemplo las siguientes funciones de forma:
[[code]]
N{3} =
                         1/2 L3 (3 L3 - 1) (3 L3 - 2)

N{7} =
                             9/2 L2 L3 (3 L3 - 1)

N{10} = 
                                  27 L1 L2 L3
[[code]]

Siendo los gráficos de estas funciones:
[[image:c5_func_forma_N3_triang_10_nodos.png]][[image:c5_func_forma_N7_triang_10_nodos.png]][[image:c5_func_forma_N10_triang_10_nodos.png]]

==Efecto del Jacobiano en la transformación isoparamétrica==
Código de MATLAB: [[file:c5_jacobiano_isoparametrico.m]] (NOTA: falta hacer los comentarios respectivos)
[[image:c5_jacobiano_isoparametrico.png]]

==Cuadratura de Gauss-Legendre para un elemento triangular==
NOTA: La siguiente información lo tomé de: [[http://math2.uncc.edu/~shaodeng/TEACHING/math5172/2010Spring/]]
[[file:c5_cuadratura_triangle.zip]]


==Programa para calcular los esfuerzos, deformaciones y desplazamientos de un sólido bidimensional utilizando elementos finitos isoparamétricos rectangulares de ocho nodos==
Calcule los campos de esfuerzos, desplazamientos y deformaciones de la estructura mostrada:
[[image:c5_isoparametric_cuad_8_nodos.jpg]]

Código de MATLAB: [[file:c5_ejemplo_isoparametricos_rect_8_nodos.zip]] 
Código de JULIA 0.5.1 (experimental): [[file:c5_ejemplo_isoparametricos_rect_8_nodos_julia_0.51.zip]]

Los esfuerzos de sigma_1 y sigma_2 calculados por el programa son: 
[[image:c5_ejemplo_isoparametricos_rect_8_nodos_s1_s2.png]]

NOTA 1: la matriz para extrapolar los los esfuerzos desde los puntos de Gauss hacia los nodos utilizando el programa [[file:c5_extrapolacion_esfuerzos.m]]
NOTA 2: falta calcular las deformaciones principales.
NOTA 3: este programa exporta los resultados a [[@http://gid.cimne.upc.es/|GiD]] por ejemplo:
[[image:c5_ejemplo_isoparametricos_rect_8_nodos_exportar_resultados_gid.png]]

Con el mismo programa, escogiendo la malla 3 (creada por Alejandro Cardona Jimenez), calculamos:
[[image:c5_gancho.png]]

==Modos de energía nula==
Con estos programas se pueden encontrar los modos de energía nula de un elemento rectangular de 4 nodos y del serendípito de 8 nodos.

[[file:c5_modos_energia_nula_rectangulo_4_nodos_K_exacta.m]] 
[[file:c5_modos_energia_nula_rectangulo_4_y_8_nodos_K_integrada.m]] (es conveniente que en este caso varíe los puntos de integración para ver los modos de energía nula que aparecen cuando se hace integración reducida)

**MODIFIQUE AMBOS PROGRAMAS PARA QUE UTILICEN EL COMANDO null() de MATLAB y así encuentre más fácilmente el espacio nulo de la matriz K**

[[image:c5_modos_energia_nula.png]]

----