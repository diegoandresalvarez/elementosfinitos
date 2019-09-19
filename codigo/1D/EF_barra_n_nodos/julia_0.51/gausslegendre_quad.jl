using Polynomials
using PyPlot

function gausslegendre_quad(m)
# Integration using a Gauss-Legendre quadrature
# Usage:
# xi, w, P = gausslegendre_quad(m)

# WHO   DATE            WHAT
# DAA   Mar 11, 2010    First algorithm in MATLAB
# DAA   Apr 04, 2017    Translation to Julia
#
# DAA - Diego Andres Alvarez Marin - daalvarez@unal.edu.co

# Gaussian quadrature xi and w

  ## Calculation of the Legendre polynomials using Bonnet's recursion:
  # (n+1) P_{n+1}(x) = (2n+1) x P_n(x) - n P_{n-1}(x)
  #         n P_n(x) = (2(n-1)+1) x P_{n-1}(x) - (n-1) P_{n-2}(x)
  #         n P_n(x) = (2n-1) x P_{n-1}(x) - (n-1) P_{n-2}(x)
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

  ## Weights: VERSION 1
  s = polyder(P[m + 1]);
  w = 2.0 ./ ((1 - xi.^2).*polyval(s, xi).^2);

  ## Weights: VERSION 2
  # A = zeros(m,m)
  # b = zeros(m,1)
  # for i=1:m
  #   A[i,:] = xi.^(i-1)'
  #   b[i]   = (1 - (-1)^i)/i
  # end
  # w = A\b;

  if ~isreal(w)
     error("m is too large. The weights cannot be complex")
  end

  return xi, w, P
end

function plot_legendre_polynomials(m = 4)
  xi, w, P = gausslegendre_quad(m)

  xx = linspace(-1,1,100)
  yy = [ polyval(P[i], xx) for i = 1:m+1 ]

  figure()
  for i = 1:m+1
    plot(xx, yy[i], label=latexstring("P_{$(i-1)}(x)"))
    plot(xi, zeros(xi), "ro")
  end
  legend()
  grid()
  axis([-1, 1, -1.1, 1.1]);
  return
end

function test_gauss_legendre_quadrature(m = 4)
  xi, w, P = gausslegendre_quad(m)

  integral(f,a,b) = ((b-a)/2)*sum(w.*f((b+a)/2 + (b-a)*xi/2))

  # OJO con este comportamiento raro de Julia
  # f1 y f2 deben tener nombres diferentes
  # si ambas funciones se llamaba igual, el resultado es erróneo
  a = 0.0; b = 0.8;
  f1(x) = 0.2 + 25x - 200x.^2 + 675x.^3 - 900x.^4 + 400x.^5
  println("Error = ", integral(f1,a,b) - 3076/1875)

  a = 0; b = pi/2
  f2(x) = sin(x)
  println("Error = ", integral(f2,a,b) - 1)
end

m = 4
plot_legendre_polynomials(m)
test_gauss_legendre_quadrature(m)
