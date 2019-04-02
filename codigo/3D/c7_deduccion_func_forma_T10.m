%%
clear, clc, close all

syms L1 L2 L3 L4

disp('Funciones de forma del tetraedro de 10 nodos')

IJKL = [ 2 0 0 0
         1 1 0 0
         0 2 0 0
         0 1 1 0
         0 0 2 0
         1 0 1 0
         0 1 0 1
         0 0 1 1
         1 0 0 1
         0 0 0 2 ];
      
coord =  IJKL(:,[2 3 4])/2;

N = cell(10,1);
for i = 1:10  
   %% Calculo de las funciones de forma   
   switch IJKL(i,1)
      case 0, lI = 1;
      case 1, lI = polyfit([0  1/2   ], [0 1  ], 1);
      case 2, lI = polyfit([0  1/2  1], [0 0 1], 2);
   end;
   switch IJKL(i,2)
      case 0, lJ = 1;
      case 1, lJ = polyfit([0  1/2   ], [0 1  ], 1);
      case 2, lJ = polyfit([0  1/2  1], [0 0 1], 2);
   end;   
   switch IJKL(i,3)
      case 0, lK = 1;
      case 1, lK = polyfit([0  1/2   ], [0 1  ], 1);
      case 2, lK = polyfit([0  1/2  1], [0 0 1], 2);
   end;      
   switch IJKL(i,4)
      case 0, lL = 1;
      case 1, lL = polyfit([0  1/2   ], [0 1  ], 1);
      case 2, lL = polyfit([0  1/2  1], [0 0 1], 2);
   end;         

   %% El round es por errores de aproximacion de la funcion polyfit
   lI = round(1000*lI)/1000;      lI = factor(poly2sym(lI,L1));
   lJ = round(1000*lJ)/1000;      lJ = factor(poly2sym(lJ,L2));
   lK = round(1000*lK)/1000;      lK = factor(poly2sym(lK,L3));
   lL = round(1000*lL)/1000;      lL = factor(poly2sym(lL,L4));   
      
   N{i} = lI*lJ*lK*lL; % = lI^i(L1) * lJ^i(L2) * lK^i(L3) * lL^i(L4)
   
   %% Se imprimen las funciones de forma   
   fprintf('\n\nN{%d} =',i); pretty(N{i});
end;

%% Se calculan e imprimen las derivadas de las funciones de forma:
fprintf('\n\n-----------------------------------------------------------');
for i = 1:10
   fprintf('\ndN{%d}/dL2 = %s',i,char(diff(N{i},L2)));
end
fprintf('\n\n-----------------------------------------------------------');
for i = 1:10
   fprintf('\ndN{%d}/dL3 = %s',i,char(diff(N{i},L3)));
end
fprintf('\n\n-----------------------------------------------------------');
for i = 1:10
   fprintf('\ndN{%d}/dL4 = %s',i,char(diff(N{i},L4)));
end
fprintf('\n\n');

return;