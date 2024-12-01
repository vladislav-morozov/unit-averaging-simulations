function obwCmap = orangeBlueMaps(n1, n2, flip)
% orangeBlueMaps Creates an orange-blue-white colormap
% 
% Args:
%   n1 (int): length of the orange-blue part of the colormap.
%   n2 (int): length of the blue-white part of the colormap.
%   flip (boolean): if true, the map is flipped upside down
%       (white-blue-orange)
%
% Retuns:
%   obwCmap (matrix): orange-blue-white or white-blue-orange colormap.
%       Expressed as a matrix of RGB triples.

% Define edge colors
blue = [80, 76, 167]/255;
orange = [236, 193, 1]/255;
white = [1, 1, 1];

% Create the colormap
% Orange-blue section
cmap_1 = [linspace(orange(1), blue(1), n1)', ...
    linspace(orange(2), blue(2), n1)', ...
    linspace(orange(3), blue(3), n1)'];
% Blue-white section
cmap_2  = [linspace(blue(1), white(1), n2)', ...
    linspace(blue(2), white(2), n2)', ...
    linspace(blue(3), white(3), n2)'];

% Combine the two sections together
obwCmap = [cmap_1; cmap_2];

% Flip the colomap if necesary
if flip
   obwCmap = flipud(obwCmap); 
end
end
