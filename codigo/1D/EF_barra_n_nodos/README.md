# Elemento finito de barra de `n` nodos

## Deducción de la matriz de rigidez `K` y el vector de fuerzas nodales equivalentes `f`
A continuación se hace la deducción de la matriz de rigidez `K` y el vector de fuerzas nodales equivalentes `f` de un elemento isoparamétrico de barra lagrangiano cuadrático (de tres nodos), suponiendo que `E`, `A` y `b` son constantes en el elemento y que los nodos están igualmente espaciados.

La salida de 
* [calculo_K_y_f.m](calculo_K_y_f.m)
* [calculo_K_y_f.py](calculo_K_y_f.m)

es:
```
K = 
    ⎡14   -16   2 ⎤
A⋅E ⎢             ⎥
───⋅⎢-16  32   -16⎥
6⋅L ⎢             ⎥
    ⎣ 2   -16  14 ⎦

f = 
    ⎡1⎤
L⋅b ⎢ ⎥
───⋅⎢4⎥
 6  ⎢ ⎥
    ⎣1⎦
```


==Comparación de varios algoritmos de interpolación implementados en MATLAB==
* Código MATLAB: [[file:c3_comparing_interpolation_algorithms.m]] 


==Funciones de forma lagrangianas para EFs de 2, 3 y 4 nodos==
* [[http://es.wikipedia.org/wiki/Interpolaci%C3%B3n_polin%C3%B3mica_de_Lagrange|Interpolación polinómica de Lagrange]]

* Código MATLAB: [[file:c3_funciones_forma_lagrangianos1D.m]] 
Al ejecutar este código obtenemos no solo la función de forma sino su dibujo, por ejemplo en el caso de las funciones lagrangianas cúbicas, tenemos:
[[code]]
Funciones de Forma Lagrangianas de CUATRO nodos:


N1 = 
    (xi - 1) (3 xi - 1) (3 xi + 1)
  - ------------------------------
                  16


N2 = 
  9 (xi - 1) (3 xi - 1) (xi + 1)
  ------------------------------
                16


N3 = 
    9 (xi - 1) (3 xi + 1) (xi + 1)
  - ------------------------------
                  16


N4 = 
  (3 xi - 1) (3 xi + 1) (xi + 1)
  ------------------------------
                16
[[code]]
[[image:c3_shape_func_lagrangianas_1D.png]]


==Como manejar los errores de redondeo con MATLAB==

A veces sucede que cuando se calculan las funciones de forma aparecen unos términos racionales extraños. Esto se debe al error de redondeo. Por ejemplo:
[[code format="matlab"]]
>> format long

>> syms x

>> N3 = poly2sym(polyfit([-1 -0.5 0 0.5 1],[0 0 1 0 0],4),x);

>> pretty(N3)

                                 3
     4         1777463176020049 x              2
  4 x  - -------------------------------- - 5 x  + 1
         10141204801825835211973625643008
[[code]]

Observe que dicho término racional tiene un valor de -1.752714012540238e-16:
[[code format="matlab"]]
>> NN3 = polyfit([-1 -0.5 0 0.5 1],[0 0 1 0 0],4)
NN3 =
   4.000000000000001  -0.000000000000000  -5.000000000000000                   0   1.000000000000000

>> NN3(2)
ans =
    -1.752714012540238e-16
[[code]]

Para eliminar dicho error de redondeo debemos obligar a MATLAB a que haga dichas cifras iguales a cero y de este modo, si se procede con el cálculo, encontramos la función de forma correcta:
[[code format="matlab"]]
>> NN3(abs(NN3) < 1e-10) = 0
NN3 =
   4.000000000000001                   0  -5.000000000000000                   0   1.000000000000000

>> N3 = poly2sym(NN3,x);

>> pretty(N3)

     4      2
  4 x  - 5 x  + 1
[[code]]



==Cuadraturas de Gauss-Legendre==
* [[http://mathworld.wolfram.com/Legendre-GaussQuadrature.html|Cuadraturas de Gauss-Legendre]]. Una tabla bonita con los pesos se encuentra [[http://de.wikipedia.org/wiki/Gau%C3%9F-Quadratur|aquí]]

Recuerde que:
[[math]]
\xi = \frac{2}{b-a}(x - \frac{a+b}{2})
[[math]]
por lo tanto
[[math]]
x = \frac{a+b}{2} + \frac{b-a}{2}\xi
[[math]]
y
[[math]]
\frac{\mathrm{d}x}{\mathrm{d}\xi} = \frac{b-a}{2}
[[math]]
por lo tanto
[[math]]
\int_a^b f(x) \mathrm{d} x = \int_{-1}^{+1} \frac{b-a}{2} f\left(\frac{a+b}{2} + \frac{b-a}{2}\xi\right) \mathrm{d} \xi
[[math]]
y utilizando las cuadraturas de Gauss-Legendre la integral anterior se puede expresar como:
[[math]]
\int_a^b f(x) \mathrm{d} x \approx \frac{b-a}{2}\sum_{i=1}^m w_i f\left(\frac{a+b}{2} + \frac{b-a}{2}\xi_i\right)
[[math]]


* Código MATLAB para generar los pesos de la cuadratura: [[file:gausslegendre_quad.m]]

* Integración de:
[[math]]
\int_0^{0.8} 0.2 + 25 x - 200 x^2 + 675x^3 - 900x^4 + 400x^5 \ \mathrm{d}x
[[math]]

Solución exacta:
[[code format="matlab"]]
syms x
int(0.2 + 25*x - 200*x^2 + 675*x^3 - 900*x^4 + 400*x^5,x,0,0.8)
[[code]]
Siendo la respuesta
[[code]]
ans =
3076/1875
[[code]]

Integración con las cuadraturas de Gauss-Legendre:
[[code format="matlab"]]
err = zeros(10,1);        % Separo la memoria
a = 0; b = 0.8;           % Limites de integracion
f = @(x) 0.2 + 25*x - 200*x.^2 +675*x.^3 - 900*x.^4 + 400*x.^5; % la funcion
solucion = 3076/1875;     % la solucion exacta
for m = 1:10;             % Vario el numero de puntos de la cuadratura
   [xi,w] = gausslegendre_quad(m);  % calculo w y c de la cuadratura
   err(m) = abs(((b-a)/2)*sum(w.*f((b+a)/2 + (b-a)*xi/2)) - solucion);
end;
figure                    % creo un lienzo
plot(err)                 % grafico el error
xlabel('Numero de puntos en la cuadratura');
ylabel('Error');
title('Cuadratura de Gauss Legendre');
grid                      % pongo la rejilla
[[code]]

Salida:
[[image:03_cuadratura_poly.png]]

**Análisis de resultados:** Recuerde que la cuadratura de Gauss-Legendre de orden n integra __EXACTAMENTE__ un polinomio de grado 2n-1 o menor. Por lo tanto con una cuadratura de grado 3 o mayor se integra exactamente este polinomio (que es de cuarto orden)

* Integración de:
[[math]]
\int_0^{\pi/2} \sin x \ \mathrm{d}x
[[math]]

Solución exacta:
[[code format="matlab"]]
syms x
int(sin(x),x,0,pi/2)
[[code]]
Siendo la respuesta
[[code]]
ans =
1
[[code]]

Integración con las cuadraturas de Gauss-Legendre:
[[code format="matlab"]]
err = zeros(10,1);
a = 0; b = pi/2;
f = @(x) sin(x);
solucion = 1;
for m = 1:10;
   [xi,w] = gausslegendre_quad(m);
   err(m) = abs(((b-a)/2)*sum(w.*f((b+a)/2 + (b-a)*xi/2)) - solucion);
end;
figure
semilogy(err); % dibujo en escala logarítmica del eje Y (para apreciar el error)
xlabel('Numero de puntos en la cuadratura');
ylabel('Error');
title('Cuadratura de Gauss Legendre');
grid
[[code]]

Salida:
[[image:03_cuadratura_sin.png]]

**Análisis de resultados:** Con la cuadratura de orden 8 o mayor se obtiene un error menor de 1e-18, y  1 - 1e-18 excede la precisión de representación de decimales del computador. Por lo tanto el error en la integración el computador lo aproxima a cero (**error de truncamiento**)

==Ejemplo de una barra sometida a carga axial constante (solución con elementos finitos de tres nodos)==
[[image:c3_ejemplo_barra.png]]
* Código MATLAB que utiliza la matriz de rigidez deducida en el curso: [[file:c3_ejemplo_barra_con_carga_axial_3_nodos_k_dado.m]] <span style="color: #ff0000;">(Este programa corre también en GNU OCTAVE)</span>

* Código MATLAB que utiliza las cuadraturas de Gauss-Legendre: [[file:c3_ejemplo_barra_con_carga_axial_3_nodos_gauss_legendre.m]] <span style="color: #ff0000;">(Para ejecutar este programa se requiere el archivo [[file:gausslegendre_quad.m]])</span>
* Código de JULIA 0.5.1. (experimental): [[file:c3_ejemplo_barra_con_carga_axial_3_nodos_gauss_legendre_julia_0.51.zip]]

Las UNICAS diferencia entre ambos programas se muestran a continuación:
[[image:c3_diferencia_entre_progs.png]]

----