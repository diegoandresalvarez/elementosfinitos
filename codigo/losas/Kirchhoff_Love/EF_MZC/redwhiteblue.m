function C = redwhiteblue(vmin, vmax, m)
% REDWHITEBLUE Colormap that linearly varies between red (+) white (=0) and blue(-)
%   REDBLUE(VMIN,VMAX,M), is an M-by-3 matrix that defines a colormap.
%   Here VMIN and VMAX defines the ranges of the colors. The colors begin
%   with a pure blue, range through shades of blue to white, where white 
%   corresponds to zero in the color scale, and then it varies through 
%   shades of red to pure red.
%   REDWHITEBLUE(VMIN,VMAX), by itself, is the same length as the current 
%   figure's colormap. If no figure exists, MATLAB creates one.
%
%   Example
%     [x,y,z] = peaks;
%     z = z + 2;   % Observe more clearly that zero corresponds to white
%     figure
%     surf(x,y,z)
%     zmin = min(z(:));
%     zmax = max(z(:));
%     colormap(redwhiteblue(zmin, zmax));
%     colorbar;
%
%             colormap(redblue)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, COLORMAP, RGBPLOT.
%   Diego Andres Alvare Marin, 4th February 2021

if nargin == 2
    m = size(get(gcf, 'colormap'), 1);
elseif nargin ~= 3
    error('REDWHITEBLUE requires 2 or 3 parameters');
end

if vmin > vmax
    error('vmin should be less than vmax');
end

M = max(abs(vmin), abs(vmax));

% From [0 0 1] to [1 1 1] to [1 0 0];
color_range = linspace(vmin,vmax, m)';
R = interp1([-M 0 M],[0 1 1], color_range);
G = interp1([-M 0 M],[0 1 0], color_range);
B = interp1([-M 0 M],[1 1 0], color_range);

C = [R G B];