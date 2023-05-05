# Taller 1: modelado de una viga utilizando elementos finitos de tensión plana

Con el objeto de contrastar la teoría aprendida y la práctica mediante el uso de un software profesional de análisis estructural, se requiere hacer el análisis de los desplazamientos, diagramas de momento flector y de fuerza cortante en una viga con diferentes alturas. Se espera que el estudiante explore, comente, discuta los conceptos aprendidos en clase, los conceptos nuevos vistos en el software y que proponga soluciones a los problemas propuestos.

Calculos realizados de forma individual con el software profesional; elaboración del trabajo escrito en parejas.

Fecha de entrega: se especificará en GOOGLE CLASSROOM. Por cada 8 horas de retraso se descontará una décima de la nota final.

## El problema propuesto
Considere la viga mostrada:

<img src="figs/viga_2023a.svg" width=800/>

Esta tiene una sección rectangular de 5 cm de ancho y $h$ de alto, tiene una luz de 2 metros y está hecha de un material con un módulo de Young $E$ = 2 GPa, un coeficiente de Poisson $\nu$ = 0.30 y una densidad de 7000 kg/m³. Se deberán analizar dos configuraciones de la viga:
* $h$=0.2 m y $q$=100 kN/m.
* $h$=2.0 m y $q$=10000 kN/m.

Tenga en cuenta que la carga distribuída triangular es tal que esta vale $q$ en $x$=1 m y 0 en $x$=2 m.

Utilizando los siguientes métodos:
* Método de los elementos finitos 2D (tensión plana) con EFs de tensión plana y los programas que se encuentran en GITHUB (MATLAB o PYTHON).
* Programa de análisis estructural que usted registró en GOOGLE CLASSROOM. NOTA: no se puede tener un software similar al de alguno de sus compañeros del curso.

## Se solicita

Se solicita calcular y graficar:

* Malla de EFs generada con el software GMSH. (OBLIGATORIO. Si no se realiza se descuentan 30 décimas)

* El desplazamiento vertical $v(x,y)$ en $y = 0$. (+5 décimas)

* La variación de los esfuerzos $\sigma_x(y)$ y $\tau_{xy}(y)$ en las secciones transversales de la viga indicadas ($x$=0.50m, 0.95m, 1.50m). (+5 décimas)

* Los gráficos con escala de colores y curvas de nivel que muestren la variación de los esfuerzos (incluyendo los esfuerzos principales con sus inclinaciones), las deformaciones y los desplazamientos. (OBLIGATORIO. Si no se realiza se descuentas 30 décimas)

* El momento flector, la fuerza cortante y la fuerza axial en las secciones transversales indicadas. Para tal fin se deben emplear las ecuaciones: (+5 décimas) 
$$V(x) = - \iint_{A(x)}   \tau_{xy}(x,y,z)\ \text{d} y \text{d} z$$

$$M(x) = - \iint_{A(x)} y \sigma_{x}(x,y,z)\ \text{d} y \text{d} z$$

$$f_\text{axial}(x) = + \iint_{A(x)} \sigma_{x}(x,y,z)\ \text{d} y \text{d} z$$

* Las reacciones en el apoyo y la distribución de dichas reacciones en la altura. (+5 décimas)

* Compare las respuestas obtenidas contra las metodologías estudiadas en el curso de análisis estructural/resistencia de materiales (solución analítica): desplazamiento vertical, momentos, fuerzas cortantes, esfuerzos, etc. (+10 décimas)

* Haga un video de máximo 15 minutos explicando como modeló dicha viga con el software profesional. Cada uno de los integrantes del grupo debe resolver individualmente este punto, usando un programa diferente al resto de compañeros del curso. (+10 décimas; sin embargo, si no se realiza se descontarán 30 décimas)

* Análisis de resultados: compare todas las soluciones obtenidas anteriormente, incluyendo la del software profesional. ¿Cuales son los porcentajes de error obtenidos? ¿Cual es el método más aproximado a la solución analítica? (+25 décimas)

## Criterios de evaluación
* Por cada punto no resuelto se descontarán 20 unidades de la nota final.

* VIDEO:
  - 0.2 Modeló adecuadamente los apoyos? la estructura? el material?
  - 0.2 Interpretó adecuadamente los gráficos resultantes? Ubicó los esfuerzos/deformaciones/desplazamientos máximos y mínimos? Hace un adecuado análisis de resultados?
  - 0.6 Exploró todas las capacidades de visualización de resultados que ofrece el software?  

* TRABAJO ESCRITO:
  - 0.5 Compara adecuadamente todos los métodos empleados.
  - 1.0 Explica qué es lo que observa en los gráficos.
  - 1.0 Hace un adecuado análisis de resultados? ¿Explica el por qué de las diferencias entre los resultados?

Recuerde explicar detalladamente como varían las cantidades en el espacio, donde están las cantidades máximas y mínimas, como se relacionan unas gráficas con otras, etc. No es solo ubicar donde están los colores, o los máximos y los mínimos, sino decir, **por qué razón se produce esa coloración**, entendiendo como la estructura está cargada, está apoyada, se deforma, etc. Se sugiere [**este (descargue archivo .PDF)**](https://github.com/diegoandresalvarez/solidos/blob/master/talleres/solidos1/ejemplo_analisis_graficos.pdf) formato para presentar los resultados. Por ejemplo con γxy: ¿qué quiere decir esta deformación? ¿cómo se está comportando en este punto la estructura dado ese valor de γxy? ¿por qué razón se produce? No es solo ubicar los máximos y los mínimos de dicha cantidad.

NOTA EXTRA:
* 5 décimas por programar algo novedoso que mejore notablemente algún aspecto del código. En este caso dichas décimas se agregarán de forma individual y no grupal.

# Otros criterios y notas
* Se solicita subir todos los archivos asociados al trabajo (.XLSX, .PDF, .MP4, .MKV, etc) directamente a GOOGLE CLASSROOM. Por favor no los empaquete en un archivo .ZIP o .RAR.

* Lo solicitado se debe subir a la plataforma GOOGLE CLASSROOM en formato PDF. El video se debe subir a GOOGLE CLASSROOM, no a YouTube u otra plataforma de videos. El video debe contener un recuadrito en el cual se vea a usted exponiendo el tema.

* Se deben entregar las presentaciones utilizadas en los videos en formato PDF.

* En ocasiones, cuando se tienen puntos de singularidad, esos valores son tan altos, que terminan opacando los colores en la estructura, mostrándolos como uno solo. En este caso, se sugiere usar una opción del software que limita los colores a mostrar a un rango. 

* Informe máximo de 25 páginas.
  * No incluir en el trabajo escrito códigos de programación, excepto pequeños bloques de máximo 10 o 15 reglones, en caso de ser necesario.
  * No es necesario escribir una introducción o un marco teórico que contenga la metodología vista en clase.

* Se sugiere aprender a manejar un programa de edición de videos. Esto les facilitará grandemente su realización.

* No los voy a penalizar en caso que ustedes obtengan desplazamientos diferentes a los que deberían dar. La experiencia ha demostrado que hay programas que simplemente no funcionan adecuadamente (aunque son pocos). Sin embargo, el estudiante debe demostrar en el video que modeló correctamente la estructura.

* Active en el software de captura de pantalla la opción para ver el ratón.

* Por mala calidad en el sonido se rebajarán 0.5 unidades. Por favor use un micrófono auxiliar (por ejemplo, un manos libres) y evite usar el micrófono del portátil para hacer el video.

* Si se sube un video de mala calidad (por ejemplo 720p de calidad o inferior) se rebajará 1.0 unidad. Mínimo 1080p. Recuerde que no tenemos limitación en el almacenamiento en GOOGLE CLASSROOM. En caso que su equipo no sea capaz de hacer videos con resolución 1080p, infórmelo previamente.

* Si se sube el video a YouTube, se descontarán 20 décimas. Los videos los debe subir directamente a GOOGLE CLASSROOM.

* Si se usa un software diferente al registrado, se descontarán 30 décimas.

* Si se modela una estructura diferente a la registrada, o si se modela usando EFs 3D, se descontarán 30 décimas.

* Si no se incluye en el video un recuadro donde se donde se vea usted hablando sobre el software se descontarán 30 décimas.

* NOTA MAXIMA FINAL = 6.0

<!---
  * VIDEO 2 (máximo 20 minutos): en este video se debe hacer una reseña crítica de las capacidades teóricas y las hipótesis fundamentales que hace el programa en cuanto al **ANALISIS DE VIGAS Y PÓRTICOS** (es decir, en cuanto a la matemática interna para el cálculo de desplazamientos, diagramas de momento flector, fuerza cortante, etc). OJO: no es mostrar como se utiliza el software, sino más mirar los manuales de referencia del mismo y mostrar que teorías, hipótesis, suposiciones, capacidades y limitaciones que tiene el programa escogido y que se utilizaron para calcular la viga. Entregar, adicionalmente, el archivo PDF utilizado en la presentación de este video. Se sugiere para la presentación tomar capturas de pantalla de los manuales de referencia del programa en cuestión. OJO: no confunda esto con la información comercial. Lo que se está solicitando está dentro de los manuales de referencia. Algunos ejemplos de buenos análisis son:
     * MIDAS GEN (análisis de vigas): https://www.youtube.com/watch?v=p06pnzg2ZPg
     * STRUSOFT FEM-DESIGN (análisis de losas): https://www.youtube.com/watch?v=xxPzgIl-mEg

Se espera que cada uno lea a fondo el manual del usuario del software. No se queden con los videos de YouTube. En el manual del usuario generalmente existe información importante sobre las hipótesis de modelado que hace cada software.

* VIDEO 2: reseña crítica de las capacidades teóricas y las hipótesis fundamentales que hace el programa en cuanto al análisis de viga (40% = 2.4)
  * 0.8 Hace un recuento de las teorías que soporta el programa, haciendo recortes del manual de referencia del mismo. Explica capacidades de cálculo y teorías que utiliza el software. 
  * 0.8 Explica hipótesis fundamentales y consejos en el modelado según se detalla en el manual del programa; hace una reseña crítica de las capacidades teóricas, las limitaciones y las hipótesis fundamentales que hace el programa en cuanto al análisis de viga
  * 0.8 Hace una reseña crítica de las ventajas/capacidades y limitaciones/suposiciones que hace el programa en cuanto al análisis de vigas y pórticos
--->