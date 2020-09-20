# Ejercicio de cálculo de pórticos y cerchas

Trabajo a presentar en grupos de máximo dos personas.

Se solicita calcular los desplazamientos/rotaciones en todos los nodos, las reacciones en los apoyos y los diagramas de momento flector, fuerza cortante y fuerza axial para las barras 3, 7, 10 y 12 de la siguiente estructura:

![](figs/torre.svg)

Cada barra tiene sección circular con radio de 10 cm, E = 200 GPa, ν = 0.30 y ρ = 2400 kg/m³ (se debe tener en cuenta el peso propio de la estructura).

En este gráfico se tiene que:
* En negro se tienen los números de los nodos.
* En rojo se tienen los números de las barras.
* Para que todos trabajemos con las mismas convenciones el nodo 1 (nodo 2) local de cada barra corresponde a aquel cuya numeración global sea la más pequeña (grande). Por ejemplo, en la barra número 12, el nodo 8 corresponde al nodo 1 local y el nodo 10 corresponde al nodo 2 local.
* Las barras en azul conforman un pórtico resistente a momentos flectores.
* Las barras en rojo (11, 12, 13 y 14) están agarradas del resto de la estructura mediante unas rótulas como la mostrada en la figura adjunta (de modo que los elementos en rojo trabajen como elementos de cerchas):

![](https://www.sikla.es/sixcms/media.php/443/thumbnails/1342_Gelenk_Joi_41_T.tif.46778.png)

Nota: si su programa en MATLAB/PYTHON no hace los diagramas solicitados, debe hacerlos a mano a partir de los datos suministrados por el programa que usted hizo, ubicando los máximos y los mínimos y los valores en los extremos.

Pasos intermedios solicitados:
* Utilizar como base los programas que figuran en [esta carpeta](../../../codigo/repaso_matricial/portico_2d).
* Deducir las fórmulas asociadas a las fuerzas nodales equivalentes que se requieran en este cálculo.
* Comparar los resultados obtenidos con un programa de cálculo estructura. NOTA un cálculo correcto no debe dar una diferencia mayor que 0.01%. Hacer un video de aproximadamente 5 minutos donde muestre los pasos más representativos del cálculo.
* Incluir las imágenes/tablas que arrojó el software.

NOTA FINAL: 
 * Quienes ya hallan visto el curso Mecánica de Sólidos 2 deben hacer el cálculo utilizando la teoría de vigas de Timoshenko-Ehrenfest. Incluir la deducción de las matrices y vectores de fuerzas nodales equivalentes.
 * Quienes no hallan visto el curso Mecánica de Sólidos 2 deben hacer el cálculo utilizando la teoría de vigas de Euler-Bernoulli (recuerde activar en el software ya sea shear area = 0 o desactivar la opción "usar deformaciones por cortante").
 
Si tienen dudas, por favor hágalas en el grupo de WhatsApp del curso, no a mi WhatsApp personal.