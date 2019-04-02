**<span style="font-size: 300%">CODIGOS DE MATLAB</span>**
<span style="color: #0000ff;
font-size: 200%;">Nota: estos códigos están hechos para que funcionen con MATLAB 2013a. No los he actualizado a versiones más nuevas de MATLAB, ya que el 2013a es el MATLAB que se tiene instalado en los computadores de la Universidad Nacional de Colombia - Sede Manizales. Por lo tanto, los programas para álgebra simbólica podrían fallar si usted utiliza versiones modernas de MATLAB.</span>

[[image:http://imgs.xkcd.com/comics/ballmer_peak.png]]
Fuente: [[http://xkcd.com/323/]]

=CAPITULO 7: sólidos tridimensionales=

==Programa para verificar las funciones de forma de los elementos tetrahédricos de cuatro nodos==
Código de MATLAB: [[file:c7_deduccion_func_forma_T4.m]]
Nota: con este programa se verifica que muchos libros de elementos finitos tienen los coeficientes de la función de forma incorrectos ya que el programa de muestra que:
[[image:c7_coeficientes_func_forma_T4.png]]

==Programa para calcular las funciones de forma de los elementos tetrahédricos de diez nodos==
Código de MATLAB: [[file:c7_deduccion_func_forma_T10.m]]

Con este código se obtuvieron las siguientes funciones de forma:
[[code]]
Funciones de forma del tetraedro de 10 nodos:

N{1} =  L1 (2 L1 - 1)
N{2} =  4 L1 L2
N{3} =  L2 (2 L2 - 1)
N{4} =  4 L2 L3
N{5} =  L3 (2 L3 - 1)
N{6} =  4 L1 L3
N{7} =  4 L4 L2
N{8} =  4 L4 L3
N{9} =  4 L4 L1
N{10} =  L4 (2 L4 - 1)
[[code]]

==Programa para calcular las funciones de forma de los elementos hexagonales de 20 nodos==
Código de MATLAB: [[file:c7_deduccion_func_forma_H20.m]]

Por ejemplo, aquí se presenta la función de forma del nodo 10: 
[[code]]
N10 = ((zeta^2 - 1)*(eta - 1)*(xi + 1))/4
[[code]]
es decir:
[[image:c7_func_forma_H20_N10.png]]

==Integración por cuadraturas de Gauss sobre dominios tetrahédricos==
Código de MATLAB: [[file:gausslegendre_quad_tetra.m]]

==Programa para calcular los esfuerzos, deformaciones y desplazamientos de un sólido tridimensional utilizando elementos finitos isoparamétricos hexahédricos de veinte nodos==
Calcule los campos de esfuerzos, desplazamientos y deformaciones de la estructura mostrada:

<span style="font-size: 200%;
color: #ff0000;">FALTA IMAGEN</span>

Código de MATLAB: [[file:c7_ejemplo_H20.zip]] 

Los esfuerzos de sigma_3 calculados por el programa y graficados en GiD son: 
[[image:c7_ejemplo_H20_esf_s3.png]]

NOTA 1: la matriz para extrapolar los los esfuerzos desde los puntos de Gauss hacia los nodos utilizando el programa [[file:c7_extrapolacion_esfuerzos_H20.m]]
NOTA 2: falta calcular las deformaciones principales.

----

