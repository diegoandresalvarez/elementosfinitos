# Placas gruesas: teoría de Mindlin

En esta carpeta y sus subcarpetas se implementan los siguientes elementos

## Integración selectiva y reducida
Se implementan los elementos finitos
* [QL9](QL9_integracion_reducida)
* [Heterosis QHET](QHET_elemento_heterosis)

Adicionalmente se muestran los modos de energía nula del elemento QL9.

## Imposición de las deformaciones angulares (assumed transverse shear strain fields)

* El programa [motivacion_ANS.m](motivacion_ANS.m) se utiliza en las diapositivas para ilustrar el origen del método

* Calculo de las matrices A*inv(P)*T para los elementos finitos QLLL, QQQQ-L, QQQQ-S y TQQL (imposición del campo vectorial de deformaciones angulares): Mirar la carpeta [APm1T_QLLL___QQQQ_S___QQQQ_L___TQQL/](APm1T_QLLL___QQQQ_S___QQQQ_L___TQQL/)

* Elemento de losa de Mindlin QQQQ-L (imposición del campo vectorial de deformaciones angulares): Mirar la carpeta
[QQQQ_L_imposicion_campo_deformaciones_cortantes/](QQQQ_L_imposicion_campo_deformaciones_cortantes/)
