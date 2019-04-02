**<span style="font-size: 300%">CODIGOS DE MATLAB</span>**

[[image:http://www.jeffpalm.com/fox/fox.jpg]]
Fuente: no fuí capaz de encontrarla... esta caricatura es de FOXTROT (http://www.foxtrot.com/). Finalmente observe que falta un \n al final del printf().

----

=CAPITULO 4 - ELEMENTOS FINITOS DE FLEXIÓN DE VIGAS=
* Programa para el cálculo de las funciones de forma del elemento de viga de Euler-Bernoulli:
** Código MATLAB: [[file:c4_func_forma_euler_bernoulli.m]] **usa el toolbox de álgebra simbolica**
** Código compatible con MATLAB 2013a: [[file:c4_func_forma_euler_bernoulli_MATLAB_2013a.m]] **usa el toolbox de álgebra simbolica**

* Programa para el cálculo de las funciones de forma del elemento de viga de Euler-Bernoulli usando polinomios interpoladores de Hermite:
** Código MATLAB: [[file:c4_hermite.m]] **usa el toolbox de álgebra simbolica**
** Código compatible con MATLAB 2013a: [[file:c4_hermite_MATLAB_2013a.m]] **usa el toolbox de álgebra simbolica**

* Programa para el cálculo de las funciones de forma del elemento de viga de Timoshenko:
** EF Lineal: [[file:c4_kf_kc_timoshenko_lineal.m]] **usa el toolbox de álgebra simbolica**
** EF Lineal (código para MATLAB 2013a): [[file:c4_kf_kc_timoshenko_lineal_MATLAB_2013a.m]] **usa el toolbox de álgebra simbolica**

** EF Cuadrático: [[file:c4_kf_kc_timoshenko_cuadratico.m]] **usa el toolbox de álgebra simbolica**
** EF Cuadrático (código para MATLAB 2013a): [[file:c4_kf_kc_timoshenko_cuadratico_MATLAB_2013a.m]] **usa el toolbox de álgebra simbolica**


* Programa para comprobar que los momentos se deben calcular en los puntos de integración de Gauss-Legendre de los elementos finitos: [[file:c4_raices_pol_Legendre.m]]
[[image:c4_nodos_legendre.png]]

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

* Interpolación acoplada (linked interpolation. Ejemplo 2.2): 
** [[file:c4_interpolacion_acoplada_ej_2_2.m]] [[file:c4_interpolacion_acoplada_ej_2_2_MATLAB_2013a.m]]
** [[file:c4_ej_4_7_o_2_3.m]] [[file:c4_ej_4_7_o_2_3_MATLAB_2013a.m]]
** [[file:c4_ej_4_8_o_2_4.m]] [[file:c4_ej_4_8_o_2_4_MATLAB_2013a.m]]
** [[file:c4_ej_4_9_o_2_5.m]] [[file:c4_ej_4_9_o_2_5_MATLAB_2013a.m]]
** [[file:c4_ej_4_10_o_2_6.m]] [[file:c4_ej_4_10_o_2_6_MATLAB_2013a.m]]
** [[file:c4_ej_4_11_o_2_7.m]] [[file:c4_ej_4_11_o_2_7_MATLAB_2013a.m]]

* Imposición de un campo de deformaciones angulares para gxz:
** [[file:c4_Bc_sustitutiva_T2.m]] 
** [[file:c4_Bc_sustitutiva_T3.m]] 
** [[file:c4_ej_4_12_o_2_8.m]] [[file:c4_ej_4_12_o_2_8_MATLAB_2013a.m]]
** [[file:c4_ej_4_13_o_2_9.m]] [[file:c4_ej_4_13_o_2_9_MATLAB_2013a.m]]

Elemento finito de Timoshenko de dos nodos calculado utilizando integración exacta:
* [[file:c4_deduccion_K_viga_exacto_Timoshenko.m]] (falta hacer la versión para MATLAB 2013a)
