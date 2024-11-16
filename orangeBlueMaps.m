function orangeBlueWhiteCmap = orangeBlueMaps(n1, n2, flip)

blue = [80, 76, 167]/255;
orange = [236, 193, 1]/255;
white = [1, 1, 1];

% % Specify the number of colors
% n_1 = 150; % You can adjust this value
% n_2  = 45;

% Create the colormap
cmap_1 = [linspace(orange(1), blue(1), n1)', ...
    linspace(orange(2), blue(2), n1)', ...
    linspace(orange(3), blue(3), n1)'];

cmap_2  = [linspace(blue(1), white(1), n2)', ...
    linspace(blue(2), white(2), n2)', ...
    linspace(blue(3), white(3), n2)'];

orangeBlueWhiteCmap = [cmap_1; cmap_2];

if flip
   orangeBlueWhiteCmap = flipud(orangeBlueWhiteCmap); 
end

end
