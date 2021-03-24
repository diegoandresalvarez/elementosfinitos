# Taller de vigas: comparación de las teorías de Euler-Bernoulli, Timoshenko-Ehrenfest y el método de los elementos finitos para tensión plana.

Con el objeto de verificar la validez de las teorías de Euler-Bernoulli, de Timoshenko-Ehrenfest, se requiere hacer el análisis de los desplazamientos, diagramas de momento flector y de fuerza cortante en una viga, utilizando dichas teorías y compararlas con las obtenidas por el método de los elementos finitos para tensión plana. Se espera que el estudiante explore, comente, discuta los conceptos aprendidos en clase, los conceptos nuevos vistos en el software y que proponga soluciones a los problemas propuestos.

Trabajo de elaboración en grupos de máximo dos integrantes. Preferiblemente, uno de los estudiantes del grupo tuvo que haber cursado Mecánica de Sólidos 2 y uno de los estudiantes tuvo que haber cursado con Diego, Elementos Finitos 1.

Fecha de entrega: se especificará en GOOGLE CLASSROOM. Por cada día de retraso se tendrán -0.3 unidades en la nota final.

## El problema propuesto
Considere la viga mostrada:

<img src="figs/viga_2021a_seccion_variable.svg"/>

Dicha viga tiene una sección rectangular y está hecha de un material con un módulo de Young *E* = 23 GPa y un coeficiente de Poisson *ν* = 0.20. Asuma el espesor *b* = 0.2 m y la altura *h* = 0.2 m, 0.5 m y 1.0 m.

Con las teorías de vigas de EB y TE se solicita calcular y graficar, para cada una de las 3 alturas *h*, los diagramas de:
* Fuerza cortante *V*
* Momento de flexión *M*
* Angulo de giro de la sección transversal *θ*
* Desplazamiento vertical *w*

Para tal fin usar los métodos:
* Matricial que resuelve las ecuaciones diferenciales de EB y TE.
* El método exacto de las funciones de discontinuidad visto en Sólidos 2.
* Programa de análisis estructural que usted registró en GOOGLE CLASSROOM (no usar como software el FTOOL).
* Programa que modele sólidos utilizando el método de los elementos finitos en tensión plana.

Con el método de los EFs para tensión plana se requiere:
* Calcular en *x* = 0.95 m, 1.05 m, 1.50 m, 1.95 m, 2.85 m, 3.00 m, 3.50 m y 3.95 m el desplazamiento horizontal *u*, esfuerzo normal *σₓ* y el esfuerzo cortante τxz. Si lo puede hacer en más puntos, aún mejor.
* A partir de *σₓ* y *τxz* estimar el momento flector *M* y la fuerza cortante *V* en esos puntos.
* Calcular en *z* = 0 m desplazamiento vertical del eje neutro de la viga *w*.

## Lo solicitado en el informe
Hacer un informe donde se:
* Comparen todas las respuestas; recuerde hacer cálculos de los porcentajes de error para mirar las diferencias entre las respuestas.
* Haga diagramas que comparen los resultados obtenidos entre ambas teorías de vigas. ¿Cuál método calculó las reacciones, momentos de flexión, fuerzas cortantes y desplazamientos más altos y más pequeños? ¿Como varían los esfuerzos *σₓ* y *τxz* al interior de la viga?
* Verifique qué tan válidas son las hipótesis de las teorías de Euler-Bernoulli y Timoshenko-Ehrenfest y la fórmula de Collingnon-Jourawski. Esto se verifica comparando contra el método de los EFs para tensión plana.
* Configure su software de modo que se emplee la misma convención para mostrar los diagramas de fuerzas cortantes y momentos flectores empleados en clase. Los momentos son positivos cuando la fibra a tracción está a compresión. El eje dependiente *M(x)* se grafica hacia arriba.

## Material a entregar
Lo solicitado se debe subir a la plataforma GOOGLE CLASSROOM en formato PDF. El video se debe subir a GOOGLE CLASSROOM, no a YouTube u otra plataforma de videos.

* VIDEO: Hacer un video de no más de 25 minutos que ilustre como resolvió el ejercicio. En el mismo video mostrar la comparación de los resultados obtenidos con MATLAB/PYTHON y con el programa escogido. No hay que hacer el análisis de resultados en el video. Esto lo hará en el trabajo escrito.
<!---
* VIDEO 2: Hacer un video de no más de 30 minutos donde se haga una reseña crítica de las capacidades teóricas y las hipótesis fundamentales que hace el programa en cuanto al análisis de vigas. OJO: no es mostrar como se utiliza el software, sino más mirar los manuales de referencia del mismo y mostrar que teorías, hipótesis, suposiciones, capacidades y limitaciones que tiene el programa escogido. Entregar, adicionalmente, el archivo PDF utilizado en la presentación de este video. En ese PDF se pueden incluir pantallazos de los manuales de referencia del software escogido. Ejemplos de excelentes videos son:
  * MIDAS GEN (análisis de vigas): https://www.youtube.com/watch?v=p06pnzg2ZPg
  * STRUSOFT FEM-DESIGN (análisis de losas): https://www.youtube.com/watch?v=xxPzgIl-mEg
--->
* Informe del trabajo con el análisis de resultados.
* Envíe, adicionalmente, los archivos de MAXIMA, EXCEL y del software empleado asociados a este ejercicio.

## Criterios de calificación
### Trabajo principal
* Calcula y grafica reacciones, *V*, *M*, *θ*, *v* con (cada punto es obligatorio, por cada punto no realizado se tendrá -1.0 unidades):
  * 0.3 EB + método matricial
  * 0.3 TE + método matricial
  * 0.3 EB + método funciones de discontinuidad
  * 0.3 TE + método funciones de discontinuidad
  * VIDEO: 0.3 EB + software que calcula vigas
  * VIDEO: 0.3 TE + software que calcula vigas
  * VIDEO: 0.6 EFs de tensión plana. Debe tener en cuenta los consejos para hacer buenas mallas, no simplemente hacer una malla supertupida. Se aconseja refinar sobre la línea de la sección transversal.
  * VIDEO: 0.6 Estima a partir de *σₓ* el momento flector *M*
  * VIDEO: 0.6 Estima a partir de *τxz* y la fuerza cortante *V*.

* INFORME: Análisis de resultados (cada punto es obligatorio, por cada punto no realizado se tendrá -2.0 unidades):
  * 1.0 compara entre sí y analiza los trés métodos de EB y los tres métodos de TE.
  * Compara contra los métodos de EB y TE los resultados obtenidos con el método de los EFs para tensión plana. Los resultados deben aparecer en el mismo gráfico que aquellos estimados por la teoría de vigas (similar a como aparece en el `main.pdf`, en la sección 9.2)
     * 0.4 *σₓ* 
     * 0.4 *τxz*
     * 0.4 *u*
     * 0.4 *w*
     * 0.4 *M* y *V*
     * 0.4 reacciones en los apoyos
  * NOTA: Hacer gráficos y/o tablas comparativos. Calcule porcentajes de error y explique los motivos por los que se tienen esas diferencias. Compare reacciones en los apoyos, diagramas de *V*, *M*, *θ* y *w* y explica el porqué de las diferencias.

<!---
* VIDEO 2: reseña crítica de las capacidades teóricas y las hipótesis fundamentales que hace el programa en cuanto al análisis de viga
  * 0.5 Explica las capacidades de cálculo y teorías que utiliza el software? Reseña crítica de las capacidades teóricas, las limitaciones y las hipótesis fundamentales que hace el programa en cuanto al análisis de viga
  * 0.5 Explica hipótesis fundamentales y consejos en el modelado según se detalla en el manual del programa? Reseña crítica de las capacidades teóricas, las limitaciones y las hipótesis fundamentales que hace el programa en cuanto al análisis de viga
  * 0.5 Explica limitaciones del programa? Reseña crítica de limitaciones que hace el programa en cuanto al análisis de viga
--->

### Otros criterios y notas
* Active en el software de captura de pantalla la opción para ver el ratón.

* Por mala calidad en el sonido se rebajarán 0.5 unidades. Por favor use un micrófono auxiliar (por ejemplo, un manos libres) y evite usar el micrófono del portátil para hacer el video.

* Si se sube un video de mala calidad (por ejemplo 720p de calidad o inferior) se rebajará 1.0 unidad. Mínimo 1080p. Recuerde que no tenemos limitación en el almacenamiento en GOOGLE CLASSROOM. En caso que su equipo no sea capaz de hacer videos con resolución 1080p, infórmelo previamente.

* Por cada día de retrazo se descontarán 3 décimas de la nota final.

* Si modela la estructura como 3D a pesar que es una de tensión/deformación plana se tendrá menos 1.0 unidad. Se debe usar necesariamente la funcionalidad de tensión/deformación plana del programa de elementos finitos (excepto si puede demostrar que este software no tiene esa opción).

* Si se sube el video a YouTube, se tendrá menos 2.0 unidades. Los videos los debe subir directamente a GOOGLE CLASSROOM.

* Si se usa un software diferente al registrado, se tendrá menos 3.0 unidades.

* Si se modela una estructura diferente a la registrada, se tendrá menos 3.0 unidades.

* Si no se incluye en el video un recuadro donde se donde se vea usted hablando sobre el software se tendrá menos 3.0 unidades.