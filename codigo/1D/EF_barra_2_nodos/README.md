**<span style="font-size: 300%">CODIGOS DE MATLAB</span>**
<span style="color: #0000ff;
font-size: 200%;">Nota: estos códigos están hechos para que funcionen con MATLAB 2013a. No los he actualizado a versiones más nuevas de MATLAB, ya que el 2013a es el MATLAB que se tiene instalado en los computadores de la Universidad Nacional de Colombia - Sede Manizales. Por lo tanto, los programas para álgebra simbólica podrían fallar si usted utiliza versiones modernas de MATLAB.</span>

[[image:http://imgs.xkcd.com/comics/ballmer_peak.png]]
Fuente: [[http://xkcd.com/323/]]

----

=CAPITULO 2=

==Cálculo ecuación 2.10 Oñate==

[[code format="matlab"]]
syms u1 u2 x x1 x2 a1 a0;
r = solve('u1 = a1*x1 + a0','u2 = a1*x2 + a0','a0','a1');
disp('a0 = '); pretty(r.a0)
disp('a1 = '); pretty(r.a1)

u = r.a1*x + r.a0;      % Se define ahora u(x) ya que conocemos a1 y a0
u = collect(u, u2);     % Se factoriza u2
u = collect(u, u1);     % Se factoriza u1
disp('u = '); pretty(u) % Observe aquí las funciones de forma
[[code]]
siendo la salida de este:
[[code]]
          u1 x2 - u2 x1
a0 =   - -------------
            x1 - x2

       u1 - u2
a1 =   -------
       x1 - x2

      /    x        x2    \      /   x1         x    \
u =   | ------- - ------- | u1 + | ------- - ------- | u2
      \ x1 - x2   x1 - x2 /      \ x1 - x2   x1 - x2 /
[[code]]


==Cálculo ecuación 2.78 Oñate==
[[code format="matlab"]]
syms x x1 x2 E A L b;          % definicion de las variables simbolicas
 
% Defino las funciones de forma
x2 = x1 + L;
N1 = (x2-x)/L;                 N2 = (x-x1)/L;

N = [N1 N2];                   % matriz de funciones de forma
B = [diff(N1,x) diff(N2,x)];   % matriz de deformación
D = E*A;                       % matriz constitutiva

% Matriz de rigidez (ecuación 2.83)
K = int(B.'*D*B, x, x1, x2);
disp('K = '); pretty(K)
 
% Vector de fuerzas nodales equivalentes (ecuación 2.83)
f = int(N.'*b, x, x1, x2);
disp('f = '); pretty(f)
[[code]]
siendo la salida de este:
[[code]]
K = 

  +-              -+
  |   E A     E A  |
  |   ---,  - ---  |
  |    L       L   |
  |                |
  |    E A   E A   |
  |  - ---,  ---   |
  |     L     L    |
  +-              -+
f = 

  +-     -+
  |  L b  |
  |  ---  |
  |   2   |
  |       |
  |  L b  |
  |  ---  |
  |   2   |
  +-     -+
[[code]]

==Ejemplo de una barra sometida a carga axial constante (solución con elementos finitos de dos nodos)==
[[image:c3_ejemplo_barra.png]]
* Código MATLAB: [[file:c2_ejemplo_barra_con_carga_axial.m]] 

==Solución del problema anterior con la función de MATLAB bvp4c==
* Código MATLAB: [[file:c2_ejemplo_barra_con_carga_axial_exacta_vs_bvp4c.m]]

----
