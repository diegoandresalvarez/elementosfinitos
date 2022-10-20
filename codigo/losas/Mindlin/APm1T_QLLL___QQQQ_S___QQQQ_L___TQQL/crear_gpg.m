function gpg = crear_gpg(var, idx)
% Crea las variables "var" asociadas a los puntos de colocación en "idx"
%
% Por ejemplo:
%
% gpg = crear_gpg('xi', [3 6 7 8 9])
%
% crea el vector:
%
% gpg = [ gxi3, gxi6, gxi7, gxi8, gxi9]

% numero de puntos de colocación
ngamma = length(idx);

% se separa la memoria
gpg = sym('g',[1 ngamma]);

% se crean las variables de acuerdo con los índices
for i = 1:ngamma
    gpg(i) = sym(['g' var num2str(idx(i))]);
end

return