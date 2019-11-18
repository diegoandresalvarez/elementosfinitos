function h = plot_EF_H20(e, xnod, LaG) 
% Grafica un EF hexaedrico de 20 nodos

% Se está graficando del mismo modo que GiD lo hace Y es el vertical

X=1; Y=2; Z=3;

h = line([xnod(LaG(e,[1 2 3 4 5 6 7 8 1 9 13 14 15 16 17 18 19 20 13]),Z); NaN; 
          xnod(LaG(e,[3 10 15]),Z); NaN; 
          xnod(LaG(e,[5 11 17]),Z); NaN; 
          xnod(LaG(e,[7 12 19]),Z)], ...
         [xnod(LaG(e,[1 2 3 4 5 6 7 8 1 9 13 14 15 16 17 18 19 20 13]),X); NaN; 
          xnod(LaG(e,[3 10 15]),X); NaN; 
          xnod(LaG(e,[5 11 17]),X); NaN; 
          xnod(LaG(e,[7 12 19]),X)], ...
         [xnod(LaG(e,[1 2 3 4 5 6 7 8 1 9 13 14 15 16 17 18 19 20 13]),Y); NaN; 
          xnod(LaG(e,[3 10 15]),Y); NaN; 
          xnod(LaG(e,[5 11 17]),Y); NaN; 
          xnod(LaG(e,[7 12 19]),Y)]);

return