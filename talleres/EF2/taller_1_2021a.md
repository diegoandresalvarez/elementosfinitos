# Taller de vigas: comparación de las teorías de Euler-Bernoulli, Timoshenko-Ehrenfest y el método de los elementos finitos para tensión plana.

Con el objeto de contrastar la teoría aprendida y la práctica mediante el uso de un software profesional de análisis estructural, se requiere hacer el análisis de los desplazamientos, diagramas de momento flector y de fuerza cortante en una viga, utilizando las teorías de Euler-Bernoulli, de Timoshenko-Ehrenfest y el método de los elementos finitos para tensión plana. Se espera que el estudiante explore, comente, discuta los conceptos aprendidos en clase, los conceptos nuevos vistos en el software y que proponga soluciones a los problemas propuestos.

Trabajo de elaboración en grupos de máximo dos integrantes. Preferiblemente, uno de los estudiantes del grupo tuvo que haber cursado Mecánica de Sólidos 2 y uno de los estudiantes tuvo que haber cursado con Diego Elementos Finitos 1.

Fecha de entrega: se especificará en GOOGLE CLASSROOM. Por cada día de retraso se tendrán -0.3 unidades en la nota final.

## El problema propuesto
Considere la viga mostrada:

<img src="figs/viga_2021a_seccion_variable.svg"/>

Dicha viga tiene una sección rectangular y está hecha de un material con un módulo de Young *E* = 20 GPa y un coeficiente de Poisson *ν* = 0.30. Asuma *h* = 0.2 m, 0.5 m y 1.0 m.

Con las teorías de vigas de EB y TE se solicita calcular y graficar los diagramas de:
* Fuerza cortante *V*
* Momento de flexión *M*
* Angulo de giro de la sección transversal *θ*
* Desplazamiento vertical *w*
Para tal fin usar los métodos:
* Matricial que resuelve la ecuación diferencial
* El método exacto de las funciones de discontinuidad visto en Sólidos 2.
* Programa de análisis estructural que usted registró en GOOGLE CLASSROOM (no usar como software el FTOOL).

Con el método de los EFs para tensión plana se requiere:
* Calcular en *x* = ***, *** y *** el desplazamiento horizontal *u*, esfuerzo normal *σₓ* y el esfuerzo cortante τxz.
* A partir de *σₓ* y τxz estimar el momento flector *M* y la fuerza cortante *V* en esos puntos.
* Calcular en *z* = 0 m desplazamiento vertical del eje neutro de la viga *w*.

## Lo solicitado en el informe
Hacer un informe donde se:
* Compare todas las las respuestas; recuerde hacer cálculos de los porcentajes de error para mirar las diferencias entre las respuestas.
* Haga diagramas que comparen los resultados obtenidos entre ambas teorías de vigas. ¿Cuál método calculó las reacciones, momentos de flexión, fuerzas cortantes y desplazamientos más altos y más pequeños? 
* Verifique qué tan válidas son las hipótesis de las teorías de Euler-Bernoulli y Timoshenko-Ehrenfest y la fórmula de Collingnon-Jourawski.
* Configure su software de modo que se emplee la misma convención para mostrar los diagramas de fuerzas cortantes y momentos flectores empleados en clase. Los momentos son positivos cuando la fibra a tracción está a compresión. El eje dependiente *M(x)* se grafica hacia arriba.

## Material a entregar
Lo solicitado se debe subir a la plataforma GOOGLE CLASSROOM en formato PDF. El video se debe subir a GOOGLE CLASSROOM, no a YouTube u otra plataforma de videos.

* VIDEO 1: Hacer un video de no más de 15 minutos que ilustre como resolvió la viga utilizando el programa seleccionado. En el mismo video mostrar la comparación de los resultados obtenidos con MATLAB/PYTHON y con el programa escogido.

* VIDEO 2: Hacer un video de no más de 30 minutos donde se haga una reseña crítica de las capacidades teóricas y las hipótesis fundamentales que hace el programa en cuanto al análisis de vigas. OJO: no es mostrar como se utiliza el software, sino más mirar los manuales de referencia del mismo y mostrar que teorías, hipótesis, suposiciones, capacidades y limitaciones que tiene el programa escogido. Entregar, adicionalmente, el archivo PDF utilizado en la presentación de este video. En ese PDF se pueden incluir pantallazos de los manuales de referencia del software escogido. Ejemplos de excelentes videos son:
  * MIDAS GEN (análisis de vigas): https://www.youtube.com/watch?v=p06pnzg2ZPg
  * STRUSOFT FEM-DESIGN (análisis de losas): https://www.youtube.com/watch?v=xxPzgIl-mEg
  
* Informe del trabajo con el análisis de resultados.
* Envíe, adicionalmente, los archivos de MAXIMA, EXCEL y del software empleado asociados a este ejercicio.

Active en el software de captura de pantalla la opción para ver el ratón.

## Criterios de calificación
### Taller principal
* Análisis y comparación de los resultados con interpolación acoplada, solución de la ecuación diferencial y el software de EFs
  * 0.3 Calculó y graficó reacciones, *V*, *M*, *θ*, *v* con Euler-Bernoulli + método matricial. Si no se realiza se resta una unidad.
  * 0.3 Calculó y graficó reacciones, *V*, *M*, *θ*, *v* con Timoshenko-Ehrenfest + método matricial. Si no se realiza se resta una unidad.
  * 0.3 Calculó y graficó reacciones, *V*, *M*, *θ*, *v* con Euler-Bernoulli + método de funciones de discontinuidad. Si no se realiza se resta una unidad.
  * 0.3 Calculó y graficó reacciones, *V*, *M*, *θ*, *v* con Timoshenko-Ehrenfest +  + método de funciones de discontinuidad. Si no se realiza se resta una unidad.
  * 0.4 VIDEO 1: Calculó y graficó reacciones, *V*, *M*, *θ*, *v* con Euler-Bernoulli + programa de análisis estructural
  * 0.4 VIDEO 1: Calculó y graficó reacciones, *V*, *M*, *θ*, *v* con Timoshenko-Ehrenfest + programa de análisis estructural
  * 0.5 INFORME: ¿Compara respuestas entre diferentes métodos de cálculo de vigas de Euler-Bernoulli? Hace gráficos y/o tablas comparativos y los explica. Calcula porcentajes de error y explica los motivos por los que se tienen esas diferencias. Compara reacciones en los apoyos, diagramas de *V*, *M*, *θ* y *v* y explica el porqué de las diferencias.
  * 0.5 INFORME: ¿Compara respuestas entre diferentes métodos de cálculo de vigas de Timoshenko? Hace gráficos y/o tablas comparativos y los explica. Calcula porcentajes de error y explica a que los motivos por los que se tienen esas diferencias. Compara reacciones en los apoyos, diagramas de *V*, *M*, *θ* y *v* y explica el porqué de las diferencias.
  * 0.5 INFORME: ¿Compara soluciones numéricas de las teorías de Euler-Bernoulli y Timoshenko? Hace gráficos y/o tablas comparativos y los explica. Calcula porcentajes de error y explica a que los motivos por los que se tienen esas diferencias. Compara reacciones en los apoyos, diagramas de *V*, *M*, *θ* y *v* y explica el porqué de las diferencias.
  * NOTA: si no se hace el video 1, se tendrán 3 unidades menos.
  * 1.0 Modela la viga usando EFs de tensión plana y se compara los resultados de desplazamiento horizontal *u*, esfuerzo normal *σₓ*, esfuerzo cortante *τxz* para al menos cuatro secciones de la viga y desplazamiento vertical del eje neutro de la viga *w*. Debe tener en cuenta los consejos para hacer buenas mallas, no simplemente hacer una malla supertupida. Adicionalmente los resultados deben aparecer en el mismo gráfico que aquellos estimados por la teoría de vigas (similar a como aparece en el main.pdf, en la sección 9.2). Si no se realiza se resta dos unidades.
  * 1.0 Compara los resultados obtenidos con las teorías de EB y TE con las que se calcularon con el método de los EFs para tensión plana.

* VIDEO 2: reseña crítica de las capacidades teóricas y las hipótesis fundamentales que hace el programa en cuanto al análisis de viga
  * 0.5 Explica las capacidades de cálculo y teorías que utiliza el software? Reseña crítica de las capacidades teóricas, las limitaciones y las hipótesis fundamentales que hace el programa en cuanto al análisis de viga
  * 0.5 Explica hipótesis fundamentales y consejos en el modelado según se detalla en el manual del programa? Reseña crítica de las capacidades teóricas, las limitaciones y las hipótesis fundamentales que hace el programa en cuanto al análisis de viga
  * 0.5 Explica limitaciones del programa? Reseña crítica de limitaciones que hace el programa en cuanto al análisis de viga

### Otros criterios
* Por mala calidad en el sonido se rebajarán 0.5 unidades. Por favor use un micrófono auxiliar (por ejemplo, un manos libres) y evite usar el micrófono del portátil para hacer el video.
* Si se sube un video de mala calidad (360p de calidad o inferior) se rebajará 1.0 unidad. Se sugiere una resolución mínima de 720p de 24 fps y preferiblemente de 1080p. Recuerde que no tenemos limitación en el almacenamiento en GOOGLE CLASSROOM.
* Por cada día de retrazo se descontarán 3 décimas de la nota final.
* Si se usa un software diferente al registrado, se tendrá menos 2.0 unidades.