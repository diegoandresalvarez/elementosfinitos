# Taller de vigas: comparación de las teorías de Euler-Bernoulli y Timoshenko

Con el objeto de contrastar la teoría aprendida y la práctica mediante el uso de un software profesional de análisis estructural, se requiere hacer el análisis de los desplazamientos, diagramas de momento flector y de fuerza cortante en una viga, utilizando las teorías de Euler-Bernoulli y de Timoshenko. Se espera que el estudiante explore, comente, discuta los conceptos aprendidos en clase, los conceptos nuevos vistos en el software y que proponga soluciones a los problemas propuestos.

Trabajo de elaboración individual

Fecha de entrega: junio 14, 2020 a las 23:59. Por cada día de retraso se tendrán -0.3 unidades en la nota final.

## El problema propuesto
Considere la viga mostrada:

<img src="figs/viga_2020a_sin_rotula.svg"/>

Dicha viga tiene una sección rectangular. En *x*=0m, la viga tiene 10 cm de ancho y 40 cm de alto; la altura varía linealmente hasta *x*=2m donde tiene una altura de 20 cm; el tramo de viga entre *x*=2m y *x*=6m tiene una 10 cm de ancho y 20 cm de alto; la viga está hecha de un material con un módulo de Young *E* = 20 GPa y un coeficiente de Poisson *ν* = 0.30. Asuma *k₁* = *k₂* = 1000 kN/m.

Se solicita calcular y graficar los diagramas de:
* Fuerza cortante *V*
* Momento de flexión *M*
* Ángulo de giro de la viga *θ*
* Desplazamiento vertical *v*

Utilizando los siguientes métodos:
<!---
* Viga de Euler-Bernoulli (solución exacta).
--->
* Viga de Euler-Bernoulli - método matricial donde la matriz de rigidez aparece al resolver la ecuación diferencial.
* Viga de Timoshenko que resulta al utilizar la interpolación acoplada que aparece en la diapositiva 69 de [04b_EF_vigas_Timoshenko.pdf](../../diapositivas/04b_EF_vigas_Timoshenko.pdf).
* Viga de Timoshenko - método matricial donde la matriz de rigidez se calcula numéricamente en cada paso usando la función `bvp4c()` o `bvp5c()` de MATLAB o su equivalente en PYTHON (+2 puntos extra).
* Programa de análisis estructural que usted registró en http://solidos2020a.shoutwiki.com/wiki/Software_para_an%C3%A1lisis_estructural_por_elementos_finitos utilizando las teorías de Euler-Bernoulli y Timoshenko. NOTA: no usar como software el FTOOL.



## Lo solicitado en el informe
Hacer un informe donde se:
* Compare todas las las respuestas; recuerde hacer cálculos de los porcentajes de error para mirar las diferencias entre las respuestas. 
* Haga diagramas que comparen los resultados obtenidos entre ambas teorías de vigas. ¿Cuál método calculó las reacciones, momentos de flexión, fuerzas cortantes y desplazamientos más altos y más pequeños? 
* Configure su software de modo que se emplee la misma convención para mostrar los diagramas de fuerzas cortantes y momentos flectores empleados en clase. Los momentos son positivos cuando la fibra a tracción está a compresión. El eje dependiente *M(x)* se grafica hacia arriba.



## Material a entregar
Lo solicitado se debe subir a la plataforma GOOGLE CLASSROOM en formato PDF. Los videos se deben subir a YouTube y se deben enlazar en GOOGLE CLASSROOM.

* Hacer un video de no más de 30 minutos donde se haga una revisión crítica de las capacidades teóricas y las hipótesis fundamentales que hace el programa en cuanto al análisis de vigas. OJO: no es mostrar como se utiliza el software, sino más mirar los manuales de referencia del mismo y mostrar que teorías, hipótesis, suposiciones, capacidades y limitaciones que tiene el programa escogido. Entregar, adicionalmente, el archivo PDF utilizado en la presentación de este video. En ese PDF se pueden incluir pantallazos de los manuales de referencia del software escogido.
* Hacer un video de no más de 15 minutos que ilustre como resolvió la viga utilizando el programa seleccionado. En el mismo video mostrar la comparación de los resultados obtenidos con MATLAB/PYTHON y con el programa escogido. 
* Informe del trabajo con el análisis de resultados.
* Envíe, adicionalmente, los archivos de MAXIMA, EXCEL y del software empleado asociados a este ejercicio.

Active en el software de captura de pantalla la opción para ver el ratón.

## Taller alternativo

* En vez de hacer este taller, se puede implementar un artículo científico sobre el tema de vigas de una revista de impacto internacional. El profesor debe dar visto bueno sobre el artículo. En este caso, la nota máxima sera 8.0. 
* Se debe realizar una exposición de máximo 20 minutos, preferiblemente hecha con [JupyterLab](https://jupyterlab.readthedocs.io/en/stable/) o [Jupyter](https://jupyter.readthedocs.io/en/latest/) o con el [MATLAB Live Editor](https://www.mathworks.com/products/matlab/live-editor.html).
* Si se modela la viga usando EFs de tensión plana y se comparan los resultados de desplazamiento horizontal *u*, esfuerzo normal *σₓ* y desplazamiento vertical del eje neutro de la viga *w*, se tendrá una unidad adicional.

Ideas: implementar viga de Reddy, Timoshenko con elementos finitos mixtos (enfoques de Hellinger-Reissner, Hu-Washizu), etc.

## Criterios de calificación
### Taller principal (nota máxima 7.0 -- es decir, podría sacar más nota, pero esta será la máxima otorgada)
* Análisis y comparación de los resultados con interpolación acoplada, solución de la ecuación diferencial y el software de EFs (60% = 3.6)*
  * 0.3 Calculó y graficó reacciones, *V*, *M*, *θ*, *v* con Euler-Bernoulli + método matricial donde la matriz de rigidez aparece al resolver la ecuación diferencial. Si no se realiza se resta una unidad.
  * 0.3 Calculó y graficó reacciones, *V*, *M*, *θ*, *v* con Timoshenko-Ehrenfest + interpolación acoplada. Si no se realiza se resta una unidad.
  * 0.3 VIDEO 1: Calculó y graficó reacciones, *V*, *M*, *θ*, *v* con Euler-Bernoulli + programa de análisis estructural
  * 0.3 VIDEO 1: Calculó y graficó reacciones, *V*, *M*, *θ*, *v* con Timoshenko-Ehrenfest + programa de análisis estructural
  * 0.8 INFORME: ¿Compara respuestas entre diferentes métodos de cálculo de vigas de Euler-Bernoulli? Hace gráficos y/o tablas comparativos y los explica. Calcula porcentajes de error y explica a que los motivos por los que se tienen esas diferencias. Compara reacciones en los apoyos, diagramas de *V*, *M*, *θ* y *v* y explica el porqué de las diferencias.
  * 0.8 INFORME: ¿Compara respuestas entre diferentes métodos de cálculo de vigas de Timoshenko? Hace gráficos y/o tablas comparativos y los explica. Calcula porcentajes de error y explica a que los motivos por los que se tienen esas diferencias. Compara reacciones en los apoyos, diagramas de *V*, *M*, *θ* y *v* y explica el porqué de las diferencias.
  * 0.8 INFORME: ¿Compara soluciones numéricas de las teorías de Euler-Bernoulli y Timoshenko? Hace gráficos y/o tablas comparativos y los explica. Calcula porcentajes de error y explica a que los motivos por los que se tienen esas diferencias. Compara reacciones en los apoyos, diagramas de *V*, *M*, *θ* y *v* y explica el porqué de las diferencias.
  * NOTA: si no se hace el video 1, se tendrán 3 unidades menos.

* VIDEO 2: revisión crítica de las capacidades teóricas y las hipótesis fundamentales que hace el programa en cuanto al análisis de viga (40% = 2.4)
  * 0.8 Explica las capacidades de cálculo y teorías que utiliza el software? Revisión crítica de las capacidades teóricas, las limitaciones y las hipótesis fundamentales que hace el programa en cuanto al análisis de viga
  * 0.8 Explica hipótesis fundamentales y consejos en el modelado según se detalla en el manual del programa? Revisión crítica de las capacidades teóricas, las limitaciones y las hipótesis fundamentales que hace el programa en cuanto al análisis de viga
  * 0.8 Explica limitaciones del programa? Revisión crítica de limitaciones que hace el programa en cuanto al análisis de viga

* EXTRAS SOLO PARA EL TALLER PRINCIPAL:
  * 1.0 Calculó y graficó reacciones, *V*, *M*, *θ*, *v* con Timoshenko-Ehrenfest + `bvp4c()` o `bvp5c()`
  * 1.0 Programa algo que mejore radicalmente los códigos de la WIKI de vigas

### Taller alterno (nota máxima 8.0 -- es decir, podría sacar más nota, pero esta será la máxima otorgada)
* 7.0 Implementa un artículo científico sobre el tema de vigas de una revista de impacto internacional. Realizar una exposición de máximo 20 minutos, preferiblemente hecha con JupyterLab o Jupyter o con el MATLAB Live Editor.

### Extra para ambos talleres
* 2.0 Modela la viga usando EFs de tensión plana y se compara los resultados de desplazamiento horizontal *u*, esfuerzo normal *σₓ*, esfuerzo cortante *τxy* para al menos cuatro secciones de la viga y desplazamiento vertical del eje neutro de la viga *w*. Debe tener en cuenta los consejos para hacer buenas mallas, no simplemente hacer una malla supertupida. Adicionalmente los resultados deben aparecer en el mismo gráfico que aquellos estimados por la teoría de vigas (similar a como aparece en el main.pdf, en la sección 9.2).

### Otros criterios
* Por mala calidad en el sonido se rebajarán 0.5 unidades. Por favor use un micrófono auxiliar (por ejemplo, un manos libres) y evite usar el micrófono del portátil para hacer el video.
* Si se sube un video de mala calidad (360p de calidad o inferior) se rebajará 1.0 unidad. Se sugiere una resolución mínima de 720p de 24 fps y preferiblemente de 1080p. Recuerde que no tenemos limitación en el almacenamiento en GOOGLE CLASSROOM.
* Por cada día de retrazo se descontarán 3 décimas de la nota final.
* Si se usa un software diferente al registrado, se tendrá menos 2.0 unidades.
