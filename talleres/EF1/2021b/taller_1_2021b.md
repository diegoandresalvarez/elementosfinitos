# Ejercicio de elementos finitos 1D: barra doblemente empotrada

NOTA MAXIMA = 5.0

Trabajo a presentar en grupos de máximo dos personas.

Informe máximo de 10 páginas. NOTA: no incluir en el trabajo escrito códigos de programación, excepto pequeños bloques de máximo 10 o 15 reglones, en caso de ser necesario.

Considere la barra doblemente empotrada, de módulo de elasticidad *E* constante, mostrada en la figura:

![](figs/barra_sec_variable.svg)

Suponga que esta barra tiene una sección transversal circular, módulo de elasticidad *E*=200 MPa, y que la carga distribuida que actúa sobre esta (no mostrada) está dada por la parábola *b*(*x*) = 5 kN/m³ *x*² para *x* ∈ [0 m, 3 m].

Se solicita calcular las fuerzas axiales, esfuerzos, deformaciones y desplazamientos en todos los puntos de la barra. Resuelva usando MATLAB o PYTHON.

1. Haga un programa para calcular la matriz de rigidez ***K*** y el vector de fuerzas nodales equivalentes ***f*** para elementos finitos con sección transversal [cónica truncada](http://es.wikipedia.org/wiki/Tronco_de_cono) y bases de radio izquierdo *r*₁ (en *x* = *x*₁) y radio derecho *r*₂ (en *x* = *x*₂). Aquí se deben deducir las fórmulas para este tipo especial de elemento finito.

2. Haga un programa para resolver el ejemplo mostrado usando 6 elementos finitos de igual longitud y la formulación deducida en el punto anterior (+0.5 unidades).

3. Vuelva a calcular, pero esta vez usando los EFs clásicos de área transversal constante (los vistos en clase) (+0.5 unidades).

4. Deduzca cual es la ecuación diferencial con sus correspondientes condiciones de frontera que describen el desplazamiento de la barra.

5. Utilizando la ecuación diferencial con sus correspondientes condiciones de frontera y las funciones `bvp4c()` o `bvp5c()` de MATLAB o la función `solve_bvp()` de PYTHON resuelva el problema (+0.5 unidades)..

6. Use un software de EFs profesional para calcular el problema anterior. Haga un video de máximo 10 minutos explicando como modeló dicha barra con el software profesional. Cada uno de los integrantes del grupo debe resolver individualmente este punto, usando un programa diferente al resto de compañeros del curso (+0.5 unidades).

7. Compare todas las soluciones obtenidas anteriormente y compare las respuestas contra la solución exacta (punto 5).

NOTAS:
* Si tienen dudas, por favor hágalas en el grupo de WhatsApp del curso, no a mi WhatsApp personal.

RÚBRICA:
* Por cada punto no resuelto se tendrá una 1.0 unidad menos.
* Si no se realiza la parte del software profesional se tendrán 2.0 unidades menos
* Por un buen análisis de resultados se obtendrán 2.0 unidades.

