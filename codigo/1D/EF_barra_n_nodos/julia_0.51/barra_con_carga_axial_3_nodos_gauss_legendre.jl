# clear
# clc
# close all        # borro la memoria, la pantalla y las figuras

## DEFINICION DEL PROBLEMA
#=
Calcule los desplazamientos y las reacciones en el empotramiento
de la viga mostrada

| b (carga distribuida de magnitud b)
|->->->->->->->->->->->->->->->->
|====*====*====*====....====*====o-> P (carga puntual P en nodo nno)
1    2    3    4          nno-1  nno
|<----longitud L de la barra---->|   el area transversal de la barra es A
=#

using Polynomials
using PyPlot

function gausslegendre_quad(m)
  # Integration using a Gauss-Legendre quadrature

  ## Calculation of the Legendre polynomials using Bonnet's recursion:
  #           P_n(x) = ((2n-1) x P_{n-1}(x) - (n-1) P_{n-2}(x))/n

  # Remember that  JULIA does not make 0-based indexing of arrays
  P = Vector{Polynomials.Poly{Float64}}(m+1)
  P[0 + 1] =     Poly([1.0])       # P_{0}(x) = 1
  P[1 + 1] = x = Poly([0.0, 1.0])  # P_{1}(x) = x
  for n = (1 + 1):(m-1 + 1)
    P[n + 1] = ((2*n - 1)*x*P[n-1 + 1] - (n-1)*P[n-2 + 1])/n
  end

  ## Roots
  xi = sort(roots(Poly(P[m + 1])));

  ## Weights
  s = polyder(P[m + 1]);
  w = 2.0 ./ ((1 - xi.^2).*polyval(s, xi).^2);

  if ~isreal(w)
     error("m is too large. The weights cannot be complex")
  end

  return xi, w
end

# -----------------------------------------------------------------
# Se usaron tres elementos isoparametricos lagrangianos cuadraticos
# -----------------------------------------------------------------

## defino las variables
nef = 3                       # numero de elementos finitos (EF)
nno = 2*nef + 1               # numero de nodos
E   = 200.0e9   # Pa          # modulo de elasticidad de la barra
A   = (0.01)^2  # m^2         # area transversal de la barra
L   = 2.0       # m           # longitud de la barra
b   = 1000.0    # N/m         # fuerza axial aplicada sobre cada EF
P   = 250.0     # N           # carga nodal al final de la barra

xnod = linspace(0,L,nno)      # posicion de los nodos

le   = fill(L/nef, nef)       # longitud de cada EF

LaG = [1 2 3                  # definicion de EFs con respecto a nodos
       3 4 5
       5 6 7]

## Parametros de la cuadratura de Gauss-Legendre
n_int_gl = 2                  # orden de la cuadratura de Gauss-Legendre

# El comando:
xi_gl, w_gl = gausslegendre_quad(n_int_gl)
# calcula las raices (xi_gl) y los pesos (w_gl) de polinomios de Legendre

# >> [x_gl,w_gl] = gausslegendre_quad(1)
# x_gl = 0;
# w_gl = 2;
# >> [x_gl,w_gl] = gausslegendre_quad(2)
# x_gl = [  -0.577350269189626;  0.577350269189626 ];
# w_gl = [   1.000000000000000;  1.000000000000000 ];
# >> [x_gl,w_gl] = gausslegendre_quad(3)
# x_gl = [  -0.774596669241483;                  0; 0.774596669241483 ];
# w_gl = [   0.555555555555556;  0.888888888888889; 0.555555555555556 ];
# >> [x_gl,w_gl] = gausslegendre_quad(4)
# x_gl = [  -0.861136311594054; -0.339981043584857; 0.339981043584856; 0.861136311594053 ];
# w_gl = [   0.347854845137453;  0.652145154862547; 0.652145154862547; 0.347854845137453 ];

## Relacion de cargas puntuales
f = zeros(nno)    # vector de fuerzas nodales equivalentes global
f[nno] = P        # relaciono la carga puntual en el nodo "nno"

## ensamblo la matriz de rigidez global y el vector de fuerzas nodales
#  equivalentes global
K = zeros(nno,nno)  # matriz de rigidez global
De = E*A            # matriz constitutiva del elemento
for e in 1:nef      # ciclo sobre todos los elementos finitos
   idx = vec(LaG[e,:])

   Je = le[e]/2     # Jacobiano del elemento ( = dx_dxi)

   # Calculo las matrices de rigidez y el vector de fuerzas nodales
   # equivalentes del elemento
   Ke = zeros(3,3)  # 3x3
   fe = zeros(3)    # 3x1
   for m in 1:n_int_gl
      # matriz de deformacion del elemento
      xi = xi_gl[m]
      Be = (1/Je)*[xi-1/2   -2*xi   xi+1/2]

      # matriz de rigidez del elemento e
      Ke += w_gl[m]*Be'*De*Be*Je

      # matriz de funciones de forma
      N = [xi.*(xi-1)/2   (1+xi).*(1-xi)   xi.*(xi+1)/2]

      # vector de fuerzas nodales equivalentes
      fe += w_gl[m]*N'*b*Je
   end

   K[idx,idx] += Ke
   f[idx]     += fe
end

## grados de libertad del desplazamiento conocidos y desconocidos
c = [ 1 ];    d = collect(2:nno)

# f = vector de fuerzas nodales equivalentes
# q = vector de fuerzas nodales de equilibrio del elemento
# a = desplazamientos

#| qd |   | Kcc Kcd || ac |   | fd |  # recuerde que qc=0 (siempre)
#|    | = |         ||    | - |    |
#| qc |   | Kdc Kdd || ad |   | fc |

## extraigo las submatrices y especifico las cantidades conocidas
Kcc = K[c,c]; Kcd = K[c,d]; fd = f[c]
Kdc = K[d,c]; Kdd = K[d,d]; fc = f[d]

# f = vector de fuerzas nodales equivalentes
# q = vector de fuerzas nodales de equilibrio del elemento
# a = desplazamientos
ac = [ 0 ];             # desplazamientos conocidos

## resuelvo el sistema de ecuaciones
ad = Kdd\(fc - Kdc*ac)       # calculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd    # calculo fuerzas de equilibrio desconocidas
a = zeros(nno);  a[c] = ac;  a[d] = ad # desplazamientos
q = zeros(nno);  q[c] = qd             # fuerzas nodales equivalentes

## se realizan unos calculos intermedios que necesitaremos mas adelante
nint = 10                # numero de puntos donde se interpolará dentro del EF
xi = linspace(-1,1,nint) # coordenadas naturales

# matriz de funciones de forma
N = [xi.*(xi-1)/2   (1+xi).*(1-xi)   xi.*(xi+1)/2]
xx    = Vector{Any}(nef) # interpol de posiciones (geometria) en el elemento
uu    = Vector{Any}(nef) # interpol desplazamientos en el elemento
axial = Vector{Any}(nef) # fuerzas axiales en el elemento
for e in 1:nef       # ciclo sobre todas los elementos finitos
   Je = le[e]/2      # Jacobiano del elemento ( = dx_dxi)
   Be = (1/Je)*[xi-1/2  -2*xi  xi+1/2] # matriz de deformacion del elemento

   # vector de desplazamientos nodales del elemento a^{(e)}
   ae = [ a[LaG[e,1]]
          a[LaG[e,2]]
          a[LaG[e,3]] ] # = a(LaG(e,:))';

   # vector de posiciones de nodos locales
   xe = [ xnod[LaG[e,1]]
          xnod[LaG[e,2]]
          xnod[LaG[e,3]] ] # = xnod(LaG(e,:))'
   xx[e] = N*xe # interpola sobre la geometría (coord naturales a geométricas)
   uu[e] = N*ae # interpola sobre los desplazamientos

   axial[e] = De*Be*ae # fuerzas axiales en elemento finito e
end

## imprimo los resultados
# format short g
println("Desplazamientos (m) = ", a)
println("Fuerzas nodales equivalentes (N) = ", f)
println("Fuerzas nodales de equilibrio (N) = ", q)

## Grafico la solucion analitica y la solucion por el MEF
## 1) grafico los desplazamientos de la barra
uexacto(x) = (-b*x.^2/2 + (P + b*L)*x)/(E*A) # solucion analitica

x = linspace(0, L, 30)             # 30 puntos unif/ distrib. entre 0 y L

figure()                           # cree un nuevo lienzo
plot(x, uexacto(x), "rx")          # grafico solucion analitica
for e = 1:nef                      # ciclo sobre todos los EFs
   plot(xx[e], uu[e], "b-")        # grafico solucion por MEF
end
title("Comparacion de la solucion analitica con el MEF para el desplazamiento")
xlabel("Eje X (m)")
ylabel("Desplazamiento (m)")
legend(["solucion exacta", "solucion por el MEF"], loc="lower right")

## 2) grafico la carga axial de la barra
Nexacta(x) = (P + b*(L-x))         # solucion analitica para la carga axial

figure()                           # cree un nuevo lienzo
plot(x, Nexacta(x), "rx")          # grafico solucion analitica
for e = 1:nef                      # ciclo sobre todos los EFs
    plot(xx[e], axial[e], "b-")    # grafico solucion por MEF
end
title("Comparacion de la solucion analitica con el MEF para la carga axial")
xlabel("Eje X (m)")
ylabel("Carga axial (N)")
legend(["solucion exacta", "solucion por el MEF"], loc="upper right")

##
return # bye, bye!
