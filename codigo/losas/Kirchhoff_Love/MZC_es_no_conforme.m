% Programa para entender la Figura 5.9 de Onate (2013)

% Defino la funcion de forma del nodo 3 para el GDL theta_x 
% el nodo de Onate 3, es nuestro nodo 4 local.
syms xi eta
N = ((eta + 1)*(xi - 1)^2*(xi + 1))/8; % Nb(2) = N(4,2)

%% Grafico esa funcion de forma junto con sus derivadas
[XI, ETA] = meshgrid(-1:0.05:1);

f = [ N diff(N,eta), diff(N,xi) ];

figure                   % Creo un lienzo
for j = 1:3
  subplot(1,3,j);        % Divido el lienzo en 3x1 dibujos
  grid on                % creo la rejilla
  hold on;               % Para que no se sobreescriban los graficos

  NN = matlabFunction(f(j), 'Vars', {'xi','eta'});

  xlabel('$\xi$', 'FontSize',30, 'interpreter','latex'); % titulo eje X
  ylabel('$\eta$','FontSize',30, 'interpreter','latex'); % titulo eje Y
  switch j % imprimo el titulo
     case 1
        title('$\overline{N}_2(\xi,\eta)$', ...
            'FontSize',30, 'interpreter','latex')
     case 2
        title('$\frac{\partial \overline{N}_2(\xi,\eta)}{\partial \eta}$', ...
            'FontSize',30, 'interpreter','latex')
     case 3
        title('$\frac{\partial \overline{N}_2(\xi,\eta)}{\partial \xi}$', ...
            'FontSize',30, 'interpreter','latex')
  end
  mesh(XI, ETA, NN(XI,ETA),'LineWidth',2); % malla de alambre
  surf(XI, ETA, NN(XI,ETA));               % superficie
  shading interp         % se interpolan los colores
  alpha 0.3              % opacidad de la superficie
  colormap winter        % mapa de colores a utilizar

  daspect([1 1 0.6]);    % similar a axis equal pero en 3D
  view(-37.5+60,30);     % vista tridimensional
end
