# Cálculo de una estructura por el método de los elementos finitos (para tensión plana).

Se requiere hacer el cálculo de las reacciones, los desplazamientos, deformaciones, esfuerzos, esfuerzos principales, esfuerzos de von Mises, cálculo del momento flector y de fuerza cortante en las secciones A y B de la estructura mostrada, utilizando el método de los elementos finitos para tensión plana. Se espera que el estudiante explore, comente, discuta los conceptos aprendidos en clase, los conceptos nuevos vistos en el software y que proponga soluciones a los problemas propuestos.

Trabajo de elaboración en grupos de máximo dos integrantes.

Fecha de entrega: se especificará en GOOGLE CLASSROOM. Por cada ocho horas de retraso se descontará una décima de la nota final.

## El problema propuesto
Considere la estructura mostrada:

<img src="figs/estructuraTP.png"/>

Dicha estructura tiene un espesor de 50 cm y está hecha de un material con un módulo de Young *E* = 23 GPa y un coeficiente de Poisson *ν* = 0.20.

Para tal fin usar los métodos:
* De los EFs para tensión plana.
* Programa profesional de análisis estructural con el método de los elementos finitos.

Nota: En las secciones A y B se requiere estimar el momento flector *M* y la fuerza cortante *V*. Para esto, busque en su software una opción de "integración sobre superficie" en su software.

## Lo solicitado en el informe
Hacer un informe donde se:
* Calcule las fuerzas en X, Y y el momento de empotramiento en el apoyo. 
* Muestren los gráficos de los desplazamientos, deformaciones, esfuerzos, esfuerzos principales, esfuerzos de von Mises en la estructura.
* Muestre la variación de los esfuerzos *σx* y *τxy* (en coordenadas locales) en el empotramiento y las secciones A y B.
* Calcule el momento flector y la fuerza cortante (en coordenadas locales) en las secciones A y B de la estructura mostrada.
* Compare las respuestas obtenidas por ambos miembros del grupo y con el método de los EFs programado en MATLAB o PYTHON.
* Haga un buen análisis de resultados. ¿Cómo se interpretan los gráficos?

NOTA: 
* Refine adecuadamente la malla de EFs y demuestre que utilizó las funcionalidades que provee el software para generar una buena malla. Una buena malla no refina innecesariamente donde no se necesita.
* Haga un estudio de convergencia de los resultados en ciertos puntos clave.

## Material a entregar
### Un informe (grupal)
Subir a la plataforma GOOGLE CLASSROOM en formato PDF. 

### Dos videos (individuales)
Los videos solicitados se deben subir a GOOGLE CLASSROOM, no a YouTube u otra plataforma de videos. Los videos debe contener un recuadrito en el cual se vea a usted exponiendo el tema.

* VIDEO 1: Hacer un video de no más de 25 minutos que ilustre como resolvió el ejercicio. En el mismo video mostrar la comparación de los resultados obtenidos con MATLAB/PYTHON y con el programa escogido. No hay que hacer el análisis de resultados en el video. Esto lo hará en el trabajo escrito.
* VIDEO 2: Hacer un video de no más de 30 minutos donde se haga una reseña crítica de las capacidades teóricas y las hipótesis fundamentales que hace el programa en cuanto al análisis por tensión plana. OJO: no es mostrar como se utiliza el software, sino más mirar los manuales de referencia del mismo y mostrar que teorías, hipótesis, suposiciones, capacidades y limitaciones que tiene el programa escogido. Entregar, adicionalmente, el archivo PDF utilizado en la presentación de este video. En ese PDF se pueden incluir pantallazos de los manuales de referencia del software escogido. Ejemplos de excelentes videos son:
  * MIDAS GEN (análisis de vigas): https://www.youtube.com/watch?v=p06pnzg2ZPg
  * STRUSOFT FEM-DESIGN (análisis de losas): https://www.youtube.com/watch?v=xxPzgIl-mEg

## Criterios de calificación
* VIDEO 2: 
  * 0.5 Explica las capacidades de cálculo y teorías que utiliza el software. Reseña crítica de las capacidades teóricas, las limitaciones y las hipótesis fundamentales que hace el programa en cuanto al análisis de viga
  
### Otros criterios y notas
* Active en el software de captura de pantalla la opción para ver el ratón.

* Por mala calidad en el sonido se rebajarán 0.5 unidades. Por favor use un micrófono auxiliar (por ejemplo, un manos libres) y evite usar el micrófono del portátil para hacer el video.

* Si se sube un video de mala calidad (por ejemplo 720p de calidad o inferior) se rebajará 1.0 unidad. Mínimo 1080p. Recuerde que no tenemos limitación en el almacenamiento en GOOGLE CLASSROOM. En caso que su equipo no sea capaz de hacer videos con resolución 1080p, infórmelo previamente.

* Por cada 6 horas de retrazo se descontará 1 décima de la nota final.

* Si modela la estructura como 3D a pesar que es una de tensión/deformación plana se tendrá menos 1.0 unidad. Se debe usar necesariamente la funcionalidad de tensión/deformación plana del programa de elementos finitos.

* Si se sube el video a YouTube, se tendrá menos 2.0 unidades. Los videos los debe subir directamente a GOOGLE CLASSROOM.

* Si se usa un software diferente al registrado, se tendrá menos 3.0 unidades.

* Si se modela una estructura diferente a la registrada, se tendrá menos 3.0 unidades.

* Si no se incluye en el video un recuadro donde se donde se vea usted hablando sobre el software se tendrá menos 3.0 unidades.