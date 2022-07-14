%% Se dibuja la malla de elementos finitos
viga_ensayo_flexion

X = 1; Y = 2;
xnod = msh.POS(:, 1:2);
LaG = msh.QUADS8(:,[1 5 2 6 3 7 4 8]);
nno = size(xnod,1);
nef = size(LaG,1);

% Calibrar la numeracion de los EFs de modo que sea igual al GMSH:
% Identifique un EF y coloque los numeros del EF de GMSH y de MATLAB
EF_gmsh = 0;
EF_matlab = 0;

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
axis equal 

borde = msh.LINES3(:,[1 3 2]);
nbordes = length(borde);
lado = zeros(nbordes,2);
for b = 1:nbordes
    for e = 1:nef
        if     isequal(borde(b,:), LaG(e,[1 2 3]))
            lado(b,:) = [e 123]; break;
        elseif isequal(borde(b,:), LaG(e,[3 4 5]))
            lado(b,:) = [e 345]; break;      
        elseif isequal(borde(b,:), LaG(e,[5 6 7]))
            lado(b,:) = [e 567]; break;
        elseif isequal(borde(b,:), LaG(e,[7 8 1]))
            lado(b,:) = [e 781]; break;
        end
    end
end

empotra = borde(:,[1 2])'
empotra(:)