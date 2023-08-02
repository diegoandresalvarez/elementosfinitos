# Observaciones generales y temario del curso "4101210 - Aplicaciones de elementos finitos 1"

## Citas para preguntas
Únicamente solicitándolas previamente, ya sea por correo electrónico o antes/después de la clase.


## Evaluaciones y trabajos
* Exámenes
   - 20% - Exámen 1 (viernes, septiembre 15 de 2023)
   - 20% - Exámen 2 (viernes, octubre 27 de 2023)
   - 20% - Exámen 3 (viernes, diciembre 1 de 2023)
* Talleres de programación en MATLAB/Python, uso de un software profesional de elementos finitos.
   - 10% - Taller 1 
   - 10% - Taller 2 
   - 20% - Taller 3 


<!---
La fecha del examen se definirá dos semanas antes de su realización.

El curso se evaluará mediante exámenes orales y talleres, así:

* **Examen oral 1:** 20%, todo el material visto en clase, diapositivas y lecturas. Tema por definir.
* **Examen oral 2:** 20%, todo el material visto en clase, diapositivas y lecturas. Tema por definir.
* **Examen oral 3:** 20%, todo el material visto en clase, diapositivas y lecturas. Tema por definir.
* **Ejercicios de programación:** 20%, se seleccionarán al azar dos ejercicios y se evaluarán. Se calificarán qué tan completos están y si incluyen todo lo solicitado.
* **Trabajo final:** 20%, uso de un software profesional de elementos finitos.
--->

En los exámenes siempre se preguntará: teoría, demostraciones, ejercicios numéricos y ejercicios de programación. <span style="color: #ff0000;">Se permite para los exámenes traer una hoja tamaño carta en la cual ustedes pueden escribir (POR UN SOLO LADO) todas las fórmulas y comandos de MATLAB que deseen. En la hoja no se pueden ni escribir programas, ni textos explicativos, ni se pueden escribir demostraciones. Dicha hoja debe ser de elaboración personal (no se pueden traer las hojas hechas por compañeros de este o semestres pasados) y debe hacerse a mano (se prohíbe explícitamente traer fotocopias/impresiones/reducciones).</span>

<!--- 
## Examenes
Los exámenes serán orales e individuales. Se realizarán siguiendo [este](https://github.com/diegoandresalvarez/solidos/blob/master/docs/protocolo_examenes_orales.md) protocolo. En ellos, más que evaluar conceptos de memoria o verificar si el estudiante entiende la matemática detrás de las ecuaciones, se evaluará la *capacidad crítica* que se tiene al momento de emplear los conceptos aprendidos.

## Criterios de calificación de los apuntes

Se pueden presentar los apuntes en un cuaderno y/o rayando directamente sobre impresiones del libro y diapositivas:

 * Los apuntes en un cuaderno se calificarán así:
   * 5.0 Apuntes completos y de buena claridad. Incluyen no solo lo enseñado en clase y en las diapositivas, sino también el contenido que el profesor asignó como lectura en los textos guía.
   * 4.0 Apuntes de buena calidad pero parcialmente completos; hay detalles que hacen falta
   * 2.5 Apuntes mediocres e incompletos: es difícil estudiar de ellos
   * 1.0 Apuntes supermalos
   * 0.0 No hizo apuntes

 * Los apuntes sobre las impresiones del libro/diapositivas en papel se calificarán así: 
   * 5.0 Hace muchas notas en el extremo de la página que complementan o ayudan a entender el texto del libro y/o de las diapositivas. Deduce fórmulas en la margen del texto. Marca los errores que encontró en el libro. Contienen las explicaciones extra que se hacen en los videos pero que no se explican en el libro.
   * 4.0 Anotaciones adicionales de buena calidad pero parcialmente completos; hay detalles que hacen falta.
   * 1.0 Se limitó a subrayar o a marcar con resaltador. Eventualmente hay notas a mano, pero son pocas. No se evidencia que estudió con juicio las hojas.
   * 0.0 No hizo apuntes o simplemente presentó un PDF resaltado.

La razón del porqué se deben hacer las notas en papel y no electrócamente es que hay estudios que demuestran que estudiar sobre papel es más efectivo que aprender sobre una pantalla. Ver por ejemplo los artículos [1](https://www.eldiario.es/consumoclaro/consumo_digital/mejor-leer-libros-impresos-electronicos_1_3220278.html) y [2](https://www.xataka.com/otros/los-estudiantes-aprenden-mucho-mas-efectivamente-de-los-libros-impresos-que-de-pantallas-aunque-ellos-creen-lo-contrario).

* Por cada día de retrazo en la entrega de los apuntes se tendrá una décima menos.
* Si los apuntes se entregan un día antes de la fecha prevista, se tendrán dos décimas adicionales.
* Si los apuntes se entregan dos días o más días antes de la fecha prevista, se tendrán cuatro décimas adicionales.
* Durante el semestre se tendrán 30 clases aproximadamente. Al final del semestre, el conjunto de todos los apuntes se dividirá en tres grupos y de cada uno de esos grupos se seleccionará al azar uno de los apuntes. Solamente se calificarán los 3 apuntes seleccionados.
--->

## Descripción de la asignatura
En este curso se enseñará la teoría de elementos finitos para la estimación de los desplazamientos, deformaciones y esfuerzos en sólidos uni-, bi-, trimensionales y axisimétricos para materiales elásticos lineales.

## Objetivos
* Aplicar las ecuaciones básicas de la mecánica de sólidos a la solución de problemas de la ingeniería mediante la utilización del método de los elementos finitos (MEF) para el análisis de estructuras en tensión y deformación planas, estructuras de revolución y estructuras tridimensionales.
* Hacer énfasis en la programación de computadores como una herramienta para obtener soluciones numéricas de problemas cuya solución analítica es extremadamente compleja.

## Metodología
El curso se desarrollará teniendo en cuenta diferentes aspectos pedagógicos como son:
* Clases presenciales: el profesor explica los conceptos relevantes en el salón de clase.
* Realización de talleres prácticos de programación que faciliten, refuercen y aplique los conocimientos adquiridos en la parte teórica cada vez que el tema lo amerite.
* Presentación y sustentación de proyectos por parte de los estudiantes.
* Trabajo dirigido fuera de clase, ya sea individual o por grupo, por parte de los estudiantes con el propósito de afianzar los conceptos aprendidos.
* Se harán prácticas con programas de elementos finitos.

## Contenido programático de elementos finitos 1
1. CONCEPTOS INICIALES
   * Nociones de cálculo variacional.
   * Principio de la energía potencial mínima y principio del trabajo virtual
   * Breve repaso de ensamblaje matricial

2. ELEMENTOS FINITOS PARA TENSIÓN AXIAL
   * Formulación débil y fuerte.
   * Funciones de forma globales y locales
   * Elementos finitos lineales de 2, 3 y *n* nodos
   * Elemento isoparamétrico unidimensional
   * Integración numérica de Gauss-Legendre en 1D

3. ELEMENTOS FINITOS PARA TENSIÓN Y DEFORMACIÓN PLANA
   * Repaso de tensión y deformación planas
   * Elemento finito triangular de deformación constante (elemento finito CST).
   * Funciones de forma triangulares y rectangulares lagrangianas y serendípitas
   * Elementos finitos isoparamétricos en 2D
   * Integración numérica de Gauss-Legendre en 2D

4. ELEMENTOS FINITOS PARA ESTRUCTURAS TRIDIMENSIONALES
   * Repaso de ecuaciones generales de elasticidad tridimensional
   * Elemento finito tetraédrico de deformación constante.
   * Funciones de forma tetraédricas y hexaédricas lagrangianas y serendípitas
   * Elementos finitos isoparamétricos en 3D
   * Integración numérica de Gauss-Legendre en 3D

5. ELEMENTOS FINITOS PARA ESTRUCTURAS DE REVOLUCIÓN
   * Problema de elasticidad en coordenadas polares.
   * Elementos finitos triangulares y rectangulares axisimétricos
   
6. CRITERIOS PARA LA GENERACIÓN DE UNA BUENA MALLA DE ELEMENTOS FINITOS
   * Reglas básicas para la creación de una malla
   * Selección del tipo de elemento finito
   * Formas aconsejables y no aconsejables de los elementos finitos, mallas estructuradas y no estructuradas
   * Mallando en el espesor, mallando huecos
   * Tipos de degeneración de los elementos finitos
   * Criterios que ayudan a verificar la calidad de la malla: relación jacobiana, relación de aspecto, errores de alisado, relación entre radios, factor de alabeo, etc.
   * Convergencia de la solución
   * Refinación de la malla de elementos finitos (métodos h y p)

## Bibliografía
<!---
Eugenio Oñate. Cálculo de estructuras por el método de elementos finitos: análisis estático lineal. Barcelona:Centro Internacional de Métodos Numéricos en Ingeniería, CIMNE 1995. 2 edición. (en la biblioteca hay 15 ejemplares: `624.171/O59c2`).

La versión en inglés (más moderna) se puede descargar así:
--->

- [Eugenio Oñate (2009). Structural Analysis with the Finite Element Method - Linear Statics, Volume 1](https://link.springer.com/book/10.1007/978-1-4020-8733-2)
- [Eugenio Oñate (2013). Structural Analysis with the Finite Element Method - Beams, Plates and Shells, Volume 2](https://link.springer.com/book/10.1007%2F978-1-4020-8743-1)

Estos libros se pueden descargar de: 
http://bases.unal.edu.co/ `->` Springer Books `->` buscar "Oñate"

## Observaciones que se quieren dejar por escrito:
### Asistencia al curso
La puerta se cerrará 10 minutos después de haber iniciado la clase (de acuerdo con el reloj del computador del salón).

### Falta a los exámenes
Siempre que usted falte a un examen, debe haber algún documento que lo exonere de dicha inasistencia. Cuando usted por algún motivo de fuerza mayor no pueda asistir al examen, usted debe avisarle al profesor con anterioridad ya sea personalmente o por correo. En esos casos en lo posible, debe demostrarlo. Por ejemplo: si le tocó viajar a su pueblo esa semana porque algo sucedió un evento familiar de trascendencia, entonces una forma de certificar que usted viajó son los tiquetes de ida y vuelta a su pueblo. Sin una excusa o una notificación previa no se repetirán los exámenes y usted tendrá como nota un cero.

### Fraude en los exámenes o trabajos
Estos se penalizarán así:

- Nota cero en el trabajo/examen en cuestión.
- Carta a la Dirección del Departamento de Ingeniería Civil reportando el suceso.
- Se pierden adicionalmente todos los privilegios que se tienen de una calificación con notas mayores a 5.0 en todas las notas obtenidas en el semestre.

### "Minuciosamente" en los exámenes
En todos los exámenes se debe relacionar con palabras las fórmulas y motivar físicamente el por qué de un procedimiento o fórmula (es decir, se debe escribir la explicación suponiendo que usted está escribiendo un libro). Si no se hace esto, se le rebajará en ese punto en particular el 50% de la nota.
