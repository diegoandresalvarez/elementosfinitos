# Taller 2: elementos finitos bidimensionales y axisimétricos

* Fecha y hora de entrega: jueves enero 16, 2020 a las 23:59.
* Presentación individual.
* Lenguajes de programación a utilizar: MATLAB (Ossa) o PYTHON (el resto del grupo).

**Notas:** 
* Por cada día de retraso en la entrega del trabajo se les descontará 0.3 unidades de la nota final.
* El trabajo es sustentable. Si no se aprueba la sustentación se obtendrá un cero en la nota de dicho trabajo.

## Criterios de calificación
* Trabajo presentado utilizando LaTeX = +10% sobre la nota final del taller.
* Codigo sin comentarios = -1.0 por ejercicio.
* Código feo/desordenado = -1.0 por ejercicio (ya que se dificulta la legibilidad del ejercicio).
* Errores en el código = -50% del punto en cuestión.
* No interpretar información dada por el programa que usted elaboró = -1.0 por ejercicio.
* No se relacionan los resultados obtenidos en el informe final = -1.0 por ejercicio.
* Hacer algo más en el código que lo dado en clase y que mejore notablemente la presentación de los resultados = +1.0 por ejercicio.
* Cuando se utiliza un software y se compare contra este los resultados obtenidos por su programa y no se explique porqué difieren los resultados, se tendrá -0.5 por ejercicio.

## Consejos/reglas para presentar el informe
* Haga una tabla de dos columnas. En la izquierda, haga el gráfico, en la derecha, su interpretación. Explique porqué el comportamiento visto en el gráfico, localice los puntos con los valores máximos y mínimos mostrados, las zonas críticas de la estructura, y cualquier otro apunte que se considere conveniente.
* No incluya en el informe el código de su programa (solo se pueden incluir fragmentos en caso extremadamente necesario). Limítese a hacer un análisis de resultados en los informes. Incluya la deducción de las ecuaciones o formulaciones que tuvo que emplear en caso que estas no se hayan discutido en clase.
* Los trabajos se deben entregar preferiblemente de forma electrónica y en formato PDF (si lo entregan impreso que sea por ambos lados de la hoja o en hojas de reciclaje para ahorrar papel). 
* Adjuntar los códigos por correo electrónico SUPERCOMENTADOS. 
* El reporte debe incluir el análisis de resultados y cualquier otra información que usted considere necesaria.

# Ejercicio 1: estructura en tensión plana, elementos T10
Dada la estructura siguiente:

![taller_T10.jpg](figs/taller_T10.jpg)

y utilizando elementos finitos triangulares de 10 nodos, haga un programa que estime y grafique:
* Desplazamientos horizontales y verticales en cada nodo (+1 unidad)
* Las deformaciones ex, ey, gxy y los esfuerzos sx, sy, txy, s1 (con su ángulo), s2 (con su ángulo), el esfuerzo cortante máximo tmáx (con sus ángulos) y el esfuerzo de von Mises en los puntos de Gauss y en los nodos de la malla de EFs (+2 punto). Si no se realiza la parte de la extrapolación de los esfuerzos a los nodos se tendrán -2 unidades.
* Las fuerzas en los apoyos (reacciones).

Los siguientes puntos son opcionales, siempre y cuando se hayan realizado los ejercicios anteriores:
 * Si se repite el ejercicio, pero esta vez creando la malla con el programa [gmsh](http://gmsh.info/), se tendrá +2 unidades. Si se muestran los resultados utilizando gmsh se tendrá +1 unidades adicionales. 
 * Si se repite el ejercicio, pero esta vez creando la malla con el programa [GiD](https://www.gidhome.com/), se tendrá +0.5 unidades. Si se muestran los resultados utilizando GiD se tendrá +0.5 unidades adicionales. Nota: se deben crear los programas para leer los archivos de la malla o para escribir los archivos con las respuestas. Le pueden preguntar a la gente que está viendo con el profesor Paredes que les indiquen como hacer el ejercicio.
 * Si se compara la respuesta con la obtenida por el programa [GetFEM++](http://getfem.org/index.html) y se muestran los resultados con [Paraview](https://www.paraview.org/), se tendrán +2 unidades.
 * Si se compara la respuesta con la obtenida por el programa de EFs que usted registró en la WIKI se obtendrá +1 unidad.

En cada uno de los cuatro casos anteriores, hacer el video respectivo que explica el proceso de creación de la malla/modelación de la estructura y subirlo a YouTube. Incluir en el informe los archivos generados. Tener en cuenta las reglas vistas en clase para la creación de una buena malla.

Asuma:
* Módulo de elasticidad = 2 GPa.
* Coeficiente de Poisson = 0.3
* Espesor = 10 cm.
* Densidad del material = 2300 kg/m³.
* Una carga distribuída constante de magnitud 100 kN/m que está actuando sobre el borde exterior (nodos 7, 14, 21, 28, 35, ..., 63, 70) y que apunta en la dirección del origen de coordenadas.
<!--- 
X(theta) = 1 + theta, Y(theta) = theta^2, teniendo en cuenta que theta varía entre 0 y pi/2
--->

# Ejercicio 2: análisis axisimétrico: esfuerzos en un pavimento flexible

Considere un pavimento flexible. Este está formado por 4 capas de diferentes rigideces y espesores, como se muestra en la figura:

![figs/pavimento_axisimetrico.png](figs/pavimento_axisimetrico.png)

Ilustración tomada de: Sam Helwany (2007). Applied soil mechanics with ABAQUS applications. John Wiley & Sons.

Asumiendo que se le aplica una carga uniformemente distribuída de 10 kPa, distribuída en un área circular de 0.5 m de radio, calcule los desplazamientos, las deformaciones y los esfuerzos en la masa del pavimento y del suelo. Elabore una gráfica que muestre como varían dichas cantidades en el eje axial (r=0). Para tal fin utilice elementos finitos rectangulares isoparamétricos de 8 nodos.

* Si se realiza el ejercicio utilizando los programas vistos en clase y creando la malla con el programa [gmsh](http://gmsh.info/), se tendrá +4 unidades. Se deben visualizar los resultados utilizando [Paraview](https://www.paraview.org/).

 * Si se compara la respuesta con la obtenida por el programa de EFs que usted registró en la WIKI se obtendrá +2 unidad.

Analice principalmente, el papel que juega la capa de asfalto y la base en la disipación de esfuerzos en el suelo (para tal fin, asuma que dichas capas no existen y observen como cambia la distribución de los esfuerzos en el suelo).

La idea de este ejercicio es que observe que la capa de asfalto y la base, al ser más rígidas, protegen las capas subyacentes (la cuales son más suaves) de los incrementos excesivos de presión que puede producir el tráfico. Dichos excesos de carga son una de las causas que pueden producir fisuras en el pavimento.

En cada uno de los casos anteriores, hacer el video respectivo que explica el proceso de creación de la malla/modelación de la estructura y subirlo a YouTube. Incluir en el informe los archivos generados. Tener en cuenta las reglas vistas en clase para la creación de una buena malla.
