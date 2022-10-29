function graficar_A_invP_T(nod, idx, A_invP_T)
% La siguiente línea corrige bug en MATLAB 2020b en LINUX que hace que no
% funcion "surf":
% En BASH invoque a MATLAB así:
% $ export MESA_LOADER_DRIVER_OVERRIDE=i965; matlab

syms xi eta

xxi  = linspace(-1, 1, 20); 
eeta = linspace(-1, 1, 20); 
[xxi, eeta] = meshgrid(xxi, eeta);

for f = 1:2
    figure
    t = tiledlayout('flow','TileSpacing','compact');
    for NNg = (A_invP_T(f,:))
        if isequal(NNg, sym(0))
            continue
        end        
        
        nexttile
        if isequal(symvar(NNg), xi)
            Ng = matlabFunction(NNg, 'Vars', xi);
            surf(xxi, eeta, Ng(xxi))    
        elseif isequal(symvar(NNg), eta)
            Ng = matlabFunction(NNg, 'Vars', eta);
            surf(xxi, eeta, Ng(eeta))
        elseif isequal(symvar(NNg), [eta, xi])
            Ng = matlabFunction(NNg, 'Vars', [xi, eta]);
            surf(xxi, eeta, Ng(xxi, eeta))
        end
        xlabel('\xi')
        ylabel('\eta')
        if f == 1
            var = '\xi';            
        else
            var = '\eta';
        end
        title(['$\gamma_{' var '} = ' latex(NNg) '$'], ...
            'interpreter', 'latex', ...
            'FontSize', 20)
        hold on 
        h1 = plot(nod(idx(1,:),1), nod(idx(1,:),2), 'ro', 'MarkerSize', 10);
        h2 = plot(nod(idx(2,:),1), nod(idx(2,:),2), 'rx', 'MarkerSize', 10);
        legend([h1 h2], 'P. Colocación xi', 'P. Colocación eta', 'Location', 'Best')
    end
end