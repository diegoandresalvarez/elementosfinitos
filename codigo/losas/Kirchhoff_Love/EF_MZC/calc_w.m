function w = calc_w(x, y, E, nu, t, a, b, p, u, v, xi, eta)
% Funcion para calcular la deformacion de la placa rectangular
% simplemente apoyada por todos los bordes que se muestra. mostrada en la Figura:
%
%    +-----------------------------------------------------------> x
%    |                                         |             |
%    |                               u         |             |
%    |                         +-----------+   |             |
%    |                               |         |             |
%    |                   +     *************   | eta         |
%    |                   |     *************   |             |
%    |                   |     *************   |             |
%    |                   |     ***carga*****   |             |
%    |                   |     ***de********   |             |
%    |                 v |  -- ***magnitud** --+             | b
%    |                   |     ***p*********                 |
%    |                   |     *************                 |
%    |                   |     *************                 |
%    |                   |     *************                 |
%    |                   +     *************      placa      |
%    |              xi               |            de         |
%    |-------------------------------+            espesor    |
%    |                                            t          |
%    |                                                       |
%    |-------------------------------------------------------+
%    |                           a
%  y V
%
% Aqui:
% x,y      coordenadas del punto
% E,nu     modulo de elasticidad, coeficiente de Poisson
% t        espesor de la placa
% a,b      ancho placa (dir x), ancho placa (dir y)
% p        intensidad de la carga
% u,v      ancho carga (dir x), ancho carga (dir y)
% xi,eta   centro de la carga
%
% La solucion teorica se encuentra en:
% Eduard Ventsel and Theodor Krauthammer (2001)
% Thin plates and shells: theory, analysis and applications.
% Marcel Dekker : New York
% paginas 53 y 54.
%
%  La deformacion teorica en la placa esta dada por:
%  \begin{equation}
%    w(x,y) = \frac{1}{\pi^4 D}\sum_{m=1}^\infty \sum_{n=1}^\infty
%    \frac{p_{mn}}{\left(\frac{m^2}{a^2} + \frac{n^2}{b^2}\right)^2}
%    \sin\left(\frac{m \pi x}{a}\right)
%    \sin\left(\frac{n \pi y}{b}\right)
%  \end{equation}
%  donde
%  \begin{equation}
%    p_{mn} = \frac{16 p}{\pi^2 m n}
%    \sin\left(\frac{m \pi \xi}{a}\right)
%    \sin\left(\frac{n \pi \eta}{b}\right)
%    \sin\left(\frac{m \pi u}{2a}\right)
%    \sin\left(\frac{n \pi v}{2b}\right)
%  \end{equation}
%
% Nota: las formulas en
% S. Timoshenko y S. Woinowsky-Krieger (1959). Theory of plates and shells
% Mc-Graw Hill : New York.
% Seccion 29: ecuacion 130 y pag 111 ecuacion a ESTAN MALAS!!!.
% La diferencia entre ambos libros es que:
% p_{mn} = \frac{16 p}{\pi^2 m n}
% lo escriben como:
% p_{mn} = \frac{16 p}{\pi^2 m n u v}

% Por:
% Diego Andres Alvarez Marin daalvarez@unal.edu.co

D = (E*t^3)/(12*(1 - nu^2)); % rigidez a la flexion de la placa
nterm = 20; % numero de terminos de la serie a tener en cuenta

s = 0;
for m = 1:nterm
   for n = 1:nterm
      pmn = (16*p/((pi^2)*m*n))*sin(m*pi*xi/a)*sin(n*pi*eta/b)*...
         sin(m*pi*u/(2*a))*sin(n*pi*v/(2*b));

      s = s + (pmn/(((m^2)/(a^2) + (n^2)/(b^2))^2))*sin(m*pi*x/a)*sin(n*pi*y/b);
   end
end

w=s/((pi^4)*D);

if x == a || y == b % para evitar los errores por redondeo
   w = 0;
end
   
return; % bye bye
