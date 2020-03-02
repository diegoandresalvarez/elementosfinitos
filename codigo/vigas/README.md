* Programa para el cálculo de los diagramas de cortante, momento, ángulo de inclinación y desplazamiento de la viga:
[[image:c4_escamilla_ej_5_5.png width="800"]]

Solución resolviendo directamente la ecuación diferencial (usando la funcion bvp5c): 
* [[file:c4_escamilla_ej_5_5_EB_eq_diff.m]]
* [[file:c4_viga_tres_apoyos_eq_diff.m]]
NOTA: este es un método que no funciona bien con MATLAB, porque requiere un tamaño de rejilla extremadamente pequeño para dar resultados muy precisos.

Solución por el método de los elementos finitos:
> Código Euler-Bernoulli: [[file:c4_ejemplo_EB.m]]
> Código Timoshenko:[[file:c4_ejemplo_T.m]] (requiere [[file:c4_ejemplo_EB.m]] para comparar)
La comparación de los diagramas de cortante, momento, ángulo de giro y deflexión vertical para una viga de longitud 19.0m y altura h=2.0m se muestra a continuación:
[[image:eb_vs_t_h_20_1.png]]
Observe que Euler-Bernoulli proporciona una solución rígida, que en este caso no es aplicable dado que h/L>0.1

