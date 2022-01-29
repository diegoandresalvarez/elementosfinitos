// Malla no estructurada

tm =  0.1;  // Tamaño de la malla normal
tmr = 0.05;  // Tamaño de la malla refinado

  
/*
                    b1                                                                                                                
                 6-------5                         ---             
                 |       |                          |              
                 |       |                          |              
                 |       |                          |              
                 |       |                          |              
                 |       |                          | h1           
                 |       |                          |              
                 |       |                          |              
                 |       |                          |              
                 |       |                          |              
               --7       4--                       ---            
           ---/  |       |  \---                    |              
       ---/      |       |      \---                | h2           
   ---/          |       |          \---            |              
8-/--------------+-------+--------------\-3        ---             
|                                         |         |              
|                                         |         |              
|                                         |         | h3           
|                                         |         |              
|                                         |         |              
1-----------------------------------------2        ---             
                                                                   
|                   b2                    |                        
|<----------------------------------------|                        
|                                         |                        

*/                                                                  


// variables asociadas a la geometría
b1 = 0.4;
b2 = 1.55;
h1 = 1.5;
h2 = 0.17;
h3 = 0.3;

// se definen los nodos
Point(1) = {          0,        0, 0, tm};
Point(2) = {         b2,        0, 0, tm};
Point(3) = {         b2,       h3, 0, tm};
Point(4) = {b2/2 + b1/2,    h2+h3, 0, tmr};
Point(5) = {b2/2 + b1/2, h1+h2+h3, 0, tm};
Point(6) = {b2/2 - b1/2, h1+h2+h3, 0, tm};
Point(7) = {b2/2 - b1/2,    h2+h3, 0, tmr};
Point(8) = {          0,       h3, 0, tm};

// se definen las líneas
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};

// se define la superficie a mallar
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Plane Surface(1) = {1};

//Physical Line("borde_inf") = {1};
//Physical Line("borde_sup") = {5};
//Physical Surface("cimentacion") = {1};

// ajustes finales de la malla		   
Recombine Surface{1};           // usar elementos cuadriláteros en vez de triángulos
Mesh.SecondOrderIncomplete = 1; // usar elementos finitos serendípitos (incompletos)
Mesh 2;                         // generar la malla 2D
SetOrder 2;                     // se definen elementos finitos de orden 2
RenumberMeshNodes;
RenumberMeshElements;

Mesh.SurfaceFaces = 0;          // ver las "caras" de los elementos finitos 2D
Mesh.Points = 1;                // ver los nodos de la malla
Mesh.NodeLabels = 1;            // ver el número de los nodos
Mesh.SurfaceLabels = 1;         // ver el número de los EFs
General.Axes = 1;               // para que aparezcan los ejes en la interfaz

// se graba la malla
Save "cimentacion.m";           // en formato MATLAB

/*
%% Se dibuja la malla de elementos finitos
cimentacion

X = 1; Y = 2;
xnod = msh.POS(:, 1:2);
LaG = msh.QUADS8(:,[1 5 2 6 3 7 4 8]);
nno = size(xnod,1);
nef = size(LaG,1);

% Calibrar la numeracion de los EFs de modo que sea igual al GMSH:
% Identifique un EF y coloque los numeros del EF de GMSH y de MATLAB
% EF_gmsh = 258;  EF_matlab = 101;

% Active esta opción si quiere ver la malla sin mirar el GMSH;
EF_gmsh = 0;  EF_matlab = 0;

figure; hold on;
cgx = zeros(1,nef); cgy = zeros(1,nef); % almacena el centro de gravedad
for e = 1:nef
   line(xnod(LaG(e,[1:8 1]),X), xnod(LaG(e,[1:8 1]),Y));
   
   % Calculo la posicion del centro de gravedad del triangulo
   cgx(e) = mean(xnod(LaG(e,[1 3 5 7]),X));
   cgy(e) = mean(xnod(LaG(e,[1 3 5 7]),Y));
   h = text(cgx(e), cgy(e), num2str(e+EF_gmsh-EF_matlab)); 
   set(h,'Color', [1 0 0]);
end
plot(xnod(:,X), xnod(:,Y), 'bo');
text(xnod(:,X), xnod(:,Y), num2str((1:nno)'));
axis equal tight 

*/
