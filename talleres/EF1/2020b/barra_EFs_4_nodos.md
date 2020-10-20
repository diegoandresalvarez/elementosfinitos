# Ejercicio de elementos finitos 1D: barra doblemente empotrada y resuelta con elementos finitos de barra de n nodos

Trabajo a presentar en grupos de máximo dos personas.

Considere la barra doblemente empotrada de sección transversal circular y de módulo de elasticidad *E* = 200 GPa, mostrada en la figura:

![barra_seccion_constante.svg](figs/barra_seccion_constante.svg)

suponga que sobre esta barra actúa una carga distribuída (no mostrada) dada por la función:

```
b(x) = sin(x) + 0.3 cos(5x) para x ∈ [0 m, 2 m] 
```
NOTA: el argumento de las funciones sin() y cos() está dado en radianes.

* Resuelva el ejemplo mostrado usando elementos finitos de barra serendípitos de 4 nodos equiespaciados. La solución se calculará utilizando *N* elementos finitos, que el usuario podrá modificar a conveniencia. La matriz de rigidez `K` y el vector de fuerzas nodales equivalentes `f` se calcularán analíticamente. Resolver implica: calcular fuerzas axiales, esfuerzos, deformaciones y desplazamientos en todos los puntos de la barra (1 unidad).

* Resuelva el punto anterior utilizando integración numérica con cuadraturas de Gauss-Legendre, en vez de la matriz `K` y el vector `f` analíticos (2 unidades).

* **PUNTO OBLIGATORIO (si no se hace, se tendrá -1 unidad)** Utilizando la ecuación diferencial con sus correspondientes condiciones de frontera y la función `bvp4c()` de MATLAB o `solve_bvp()` de PYTHON, calcular la solución exacta (fuerzas axiales, esfuerzos, deformaciones y desplazamientos en todos los puntos de la barra) y compararla con las soluciones estimadas por el método de los EFs. ¿Cual es el error del método de los EFs? 

* Verifique la precisión del cálculo de las fuerzas axiales estimados con el punto anterior y por el método de los elementos finitos en los puntos de integración de Gauss-Legendre (o puntos de esfuerzo de Barlow). ¿Qué porcentaje de error hay en esa estimación? Haga la interpolación de las fuerzas axiales a partir de las estimaciones en dichos puntos (2 unidades).

Si tienen dudas, por favor hágalas en el grupo de WhatsApp del curso, no a mi WhatsApp personal.
